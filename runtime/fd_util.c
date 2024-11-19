#include <stdio.h>
#include <string.h>
#include <fcntl.h>
#include <dirent.h>
#include <unistd.h>
#include <limits.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <linux/kcmp.h>
#include <sys/syscall.h>
#include <unistd.h>

struct open_fd {
    int fd;
    int dup_of;
    mode_t mode;
    off_t offset;     /* only exists for mode == O_RDONLY */
    char *filename;
};

void print_open_fd(struct open_fd *fd_entry) 
{
    printf("%s %c %ld\n", fd_entry->filename,
	   fd_entry->mode == O_RDONLY ? 'r' : 'w',
	   fd_entry->offset);
}

void free_open_fd(struct open_fd *fd_entry)
{
    if (fd_entry) {
	if (fd_entry->filename) {
	    free(fd_entry->filename);
	}
        free(fd_entry);
    }
}


struct open_fd_vec {
    int len;
    int size;
    struct open_fd **open_fds;
};

struct open_fd_vec *create_open_fd_vec(void)
{
    struct open_fd_vec *vec;
    vec = (struct open_fd_vec*) malloc(sizeof(*vec));
    vec->len = 16;
    vec->size = 0;
    vec->open_fds = (struct open_fd**) malloc(vec->len * sizeof(*vec->open_fds));
    return vec;
}

int are_same_fd(int fd1, int fd2)
{
    int ret = syscall(SYS_kcmp, getpid(), getpid(), KCMP_FILE, fd1, fd2);
    return !ret;
}

void push_to_open_fd_vec_dedup(struct open_fd_vec *vec, struct open_fd *fd)
{
    int i, j;
    if (vec->size == vec->len) {
	vec->len = vec->len * 2;
	vec->open_fds = realloc(vec->open_fds, vec->len * sizeof(*vec->open_fds));
    }
    vec->open_fds[vec->size] = fd;
    i = vec->size;
    for (j = 1; j < vec->size; j++) {
	if (are_same_fd(vec->open_fds[j]->fd, vec->open_fds[i]->fd)) {
	    vec->open_fds[i]->dup_of = vec->open_fds[j]->fd;
	}
    }
    vec->size++;
}

void push_to_open_fd_vec(struct open_fd_vec *vec, struct open_fd *fd)
{
    if (vec->size == vec->len) {
	vec->len = vec->len * 2;
	vec->open_fds = realloc(vec->open_fds, vec->len * sizeof(*vec->open_fds));
    }
    vec->open_fds[vec->size] = fd;
    vec->size++;
}

void release_open_fd_vec(struct open_fd_vec *vec)
{
    int i;
    for (i = 0; i < vec->size; i++) free_open_fd(vec->open_fds[i]);
    free(vec->open_fds);
    free(vec);
}


void print_open_fd_command(struct open_fd_vec *vec)
{
    int i, j;
    for (i = 0; i < vec->size; i++) {
	if (vec->open_fds[i]->dup_of == -1) {
	    struct open_fd *open_fd = vec->open_fds[i];
	    /* no duplicate fd before */
	    printf("exec %d%c%s\n", open_fd->fd, open_fd->mode == O_RDONLY ? '<' : '>', open_fd->filename);
	} else {
	    struct open_fd *open_fd = vec->open_fds[i];
	    j = vec->open_fds[i]->dup_of;
	    printf("exec %d%c&%d\n", open_fd->fd, open_fd->mode == O_RDONLY ? '<' : '>', j);
	}
    }
}


#define BUFLEN 1024
struct open_fd *load_open_fd(int fd)
{
    char fd_path[PATH_MAX];
    char file_path[PATH_MAX + 1];
    char buf[BUFLEN];
    ssize_t len;
    struct stat statbuf;
    struct open_fd *open_fd_entry;

    open_fd_entry = (struct open_fd *)malloc(sizeof(*open_fd_entry));
    memset(open_fd_entry, 0, sizeof(*open_fd_entry));
    open_fd_entry->dup_of = -1;

    open_fd_entry->fd = fd;
    snprintf(fd_path, sizeof(fd_path), "/proc/%d/fd/%d", getpid(), fd);
    len = readlink(fd_path, file_path, sizeof(file_path) - 1);

    if (len == -1) {
	perror("readlink");
	goto err;
    }

    file_path[len] = '\0';
    open_fd_entry->filename = malloc(len + 1);
    memcpy(open_fd_entry->filename, file_path, len + 1);

    if (lstat(fd_path, &statbuf) < 0) {
	perror("fstat");
	goto err;
    }
    switch (statbuf.st_mode & (S_IWUSR | S_IRUSR)) {
    case S_IRUSR | S_IWUSR:
	switch (fd) {
        case 0:
	    open_fd_entry->mode = O_RDONLY;
	    break;
        case 1:
	case 2:
        default:
	    open_fd_entry->mode = O_WRONLY;
            break;
	    /* perror("permission"); */
	    /* goto err; */
        }
	break;
    case S_IWUSR:
	open_fd_entry->mode = O_WRONLY;
	break;
    case S_IRUSR:
	open_fd_entry->mode = O_RDONLY;
        break;
    default:
	perror("mode");
	goto err;
    }
    if (1) {
	int fdinfo_fd, i;
	char *offstr_start = NULL;
	snprintf(fd_path, sizeof(fd_path), "/proc/%d/fdinfo/%d", getpid(), fd);
	fdinfo_fd = open(fd_path, O_RDONLY);
	if (fdinfo_fd < 0) {
	    perror("open fdinfo");
	    goto err;
	}
        if (read(fdinfo_fd, buf, BUFLEN) < 0) {
	    close(fdinfo_fd);
	    goto err;
        }
        for (i = 0; i < BUFLEN; i++) {
	    if (buf[i] >= '0' && buf[i] <= '9' && offstr_start == NULL)
		offstr_start = &buf[i];
	    if (buf[i] == '\n') {
		buf[i] = '\0';
		break;
	    }
        }
	open_fd_entry->offset = strtol(offstr_start, NULL, 0);
	close(fdinfo_fd);
    }
    return open_fd_entry;
err:
    free_open_fd(open_fd_entry);
    return NULL;
}

void dump_open_fd_vec(const struct open_fd_vec *vec, const char *outfile)
{
    int i = 0, outf;
    outf = open(outfile, O_WRONLY | O_CREAT | O_TRUNC, S_IRUSR | S_IWUSR);
    if (outf < 0) {
	perror("open");
        exit(1);
    }
    for (i = 0; i < vec->size; i++) {
	struct open_fd *open_fd = vec->open_fds[i];
	if (open_fd->dup_of == -1) {
	    dprintf(outf, "%d %c %ld %s\n", open_fd->fd, open_fd->mode == O_RDONLY ? 'r' : 'w',
		   open_fd->offset, open_fd->filename);
	} else {
	    dprintf(outf, "%d %c %ld %d\n", open_fd->fd, 'd', open_fd->offset, open_fd->dup_of);
	}
    }
    close(outf);
}

struct open_fd_vec *load_open_fd_vec(const char *infile)
{
    struct open_fd_vec *vec;
    FILE *f;
    char buf[BUFLEN];
    f = fopen(infile, "r");
    if (f == NULL) {
	perror("fopen");
	exit(1);
    }
    vec = create_open_fd_vec();
    while (fgets(buf, BUFLEN, f) != NULL) {
	char mode;
	int pos;
	struct open_fd *open_fd = (struct open_fd *)malloc(sizeof(struct open_fd));
	memset(open_fd, 0, sizeof(struct open_fd));
	sscanf(buf, "%d %c %ld %n", &open_fd->fd, &mode, &open_fd->offset, &pos);
	if (mode == 'd') {
	    sscanf(&buf[pos], "%d", &open_fd->dup_of);
	} else {
	    int len = strlen(&buf[pos]);
	    open_fd->mode = (mode == 'r' ? O_RDONLY : O_WRONLY);
	    open_fd->dup_of = -1;
	    open_fd->filename = malloc(len+1);
	    strcpy(open_fd->filename, &buf[pos]);
	    open_fd->filename[len-1] = '\0';
	}
	push_to_open_fd_vec(vec, open_fd);
    }
    fclose(f);
    return vec;
}

void setup_open_fds(struct open_fd_vec *vec, const char *partial_restore_dir)
{
    int i;
    for (i = 0; i < vec->size; i++) {
	struct open_fd *open_fd = vec->open_fds[i];
	int target_fd, fd;
	target_fd = open_fd->fd;
	if (open_fd->dup_of == -1) {
	    if (open_fd->mode == O_RDONLY) {
		fd = open(open_fd->filename, open_fd->mode);
	    } else {
		char filename[PATH_MAX+16];
		sprintf(filename, "%s/%d", partial_restore_dir, target_fd);
		fd = open(filename, open_fd->mode | O_CREAT | O_TRUNC, 0644);
	    }
	    if (fd != target_fd) {
		dup2(fd, target_fd);
		close(fd);
	    }
	} else {
	    dup2(open_fd->dup_of, target_fd);
	}
    }
}

void seek_fds(const char *infile)
{
    int i;
    struct open_fd_vec *vec = load_open_fd_vec(infile);
    for (i = 0; i < vec->size; i++) {
	struct open_fd *open_fd = vec->open_fds[i];
	if (open_fd->dup_of == -1 && open_fd->mode == O_RDONLY) {
	    lseek(open_fd->fd, open_fd->offset, SEEK_SET);
	}
    }
    release_open_fd_vec(vec);
}

void save_open_file_descriptors(const char *outfile) {
    char fd_path[PATH_MAX];
    struct dirent *entry;
    DIR *dir;
    struct open_fd_vec *open_fd_vec = create_open_fd_vec();

    snprintf(fd_path, sizeof(fd_path), "/proc/%d/fd", getpid());
    dir = opendir(fd_path);

    if (dir == NULL) {
        perror("opendir");
        exit(EXIT_FAILURE);
    }

    while ((entry = readdir(dir)) != NULL) {
	struct open_fd *open_fd_entry;
        if (entry->d_name[0] == '.' || atoi(entry->d_name) == dirfd(dir))
            continue;

	open_fd_entry = load_open_fd(atoi(entry->d_name));
	push_to_open_fd_vec_dedup(open_fd_vec, open_fd_entry);
    }

    dump_open_fd_vec(open_fd_vec, outfile);
    release_open_fd_vec(open_fd_vec);
    closedir(dir);
}

void partial_restore_execve(char *filename, char *partial_restore_dir, char *argv[]) {
    char *program = argv[0];
    struct open_fd_vec *vec = load_open_fd_vec(filename);
    setup_open_fds(vec, partial_restore_dir);
    if (execvp(program, &argv[0]) == -1) {
	perror("execve");
	exit(1);
    }
}

void print_usage(const char *prog_name)
{
    printf("Usage %s [-f dumpfile] [-s|-r command| -p command]\n", prog_name);
}

#define SAVE_MODE 0
#define RESTORE_MODE 1
#define PARTIAL_RESTORE_MODE 2
#define SEEK_MODE 3

int main(int argc, char *argv[]) {
    int opt, mode = SAVE_MODE;
    char *filename = NULL;
    char partial_restore_dir[PATH_MAX];
    while ((opt = getopt(argc, argv, "srkp:f:")) != -1) {
        switch (opt) {
	case 's':
	    mode = SAVE_MODE;
	    break;
	case 'r':
	    mode = RESTORE_MODE;
	    break;
	case 'k':
	    mode = SEEK_MODE;
	    break;
	case 'p':
	    mode = PARTIAL_RESTORE_MODE;
	    if (realpath(optarg, partial_restore_dir) == NULL) {
		perror("realpath");
		return 1;
	    }
	    goto done_parse;
	case 'f':
	    filename = optarg;
	    break;
	default:
	    print_usage(argv[0]);
	    exit(EXIT_FAILURE);
        }
    }
done_parse:
    if (filename == NULL || optind == argc && (mode == RESTORE_MODE || mode == PARTIAL_RESTORE_MODE)) {
	print_usage(argv[0]);
	exit(1);
    }
    switch (mode) {
    case SAVE_MODE:
	save_open_file_descriptors(filename);
	break;
    case RESTORE_MODE:
	printf("not implemented\n");
	break;
    case PARTIAL_RESTORE_MODE:
	partial_restore_execve(filename, partial_restore_dir, &argv[optind]);
	break;
    case SEEK_MODE:
	seek_fds(filename);
	break;
    default:
	exit(1);
    }
    return 0;
}
