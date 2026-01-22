/*-
 * Copyright (c) 1991, 1993
 *	The Regents of the University of California.  All rights reserved.
 * Copyright (c) 1997-2005
 *	Herbert Xu <herbert@gondor.apana.org.au>.  All rights reserved.
 *
 * This code is derived from software contributed to Berkeley by
 * Kenneth Almquist.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 3. Neither the name of the University nor the names of its contributors
 *    may be used to endorse or promote products derived from this software
 *    without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE REGENTS AND CONTRIBUTORS ``AS IS'' AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
 * OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 */

#include <stdio.h>	/* defines BUFSIZ */
#include <fcntl.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>

/*
 * This file implements the input routines used by the parser.
 */

#include "eval.h"
#include "shell.h"
#include "redir.h"
#include "syntax.h"
#include "input.h"
#include "output.h"
#include "options.h"
#include "memalloc.h"
#include "error.h"
#include "alias.h"
#include "parser.h"
#include "main.h"
#ifndef SMALL
#include "myhistedit.h"
#endif

#define EOF_NLEFT -99		/* value of parsenleft when EOF pushed back */
#define IBUFSIZ (BUFSIZ + 1)


MKINIT struct parsefile basepf;	/* top level input file */
MKINIT char basebuf[IBUFSIZ];	/* buffer for top level input file */
struct parsefile *parsefile = &basepf;	/* current input file */
int whichprompt;		/* 1 == PS1, 2 == PS2 */

#ifndef SMALL
EditLine *el;			/* cookie for editline package */
#endif

STATIC void pushfile(void);
static int preadfd(void);
/*
static void setinputfd(int fd, int push); 
*/ // libdash
static int preadbuffer(void);

#ifdef mkinit
INCLUDE <stdio.h>
INCLUDE "input.h"
INCLUDE "error.h"

INIT {
	basepf.nextc = basepf.buf = basebuf;
	basepf.linno = 1;
}

RESET {
	/* clear input buffer */
	basepf.lleft = basepf.nleft = 0;
	popallfiles();
}
#endif


/*
 * Read a character from the script, returning PEOF on end of file.
 * Nul characters in the input are silently discarded.
 */

int
pgetc(void)
{
	int c;

	if (parsefile->unget)
		return parsefile->lastc[--parsefile->unget];

	if (--parsefile->nleft >= 0)
		c = (signed char)*parsefile->nextc++;
	else
		c = preadbuffer();

	parsefile->lastc[1] = parsefile->lastc[0];
	parsefile->lastc[0] = c;

	return c;
}


/*
 * Same as pgetc(), but ignores PEOA.
 */

int
pgetc2()
{
	int c;
	do {
		c = pgetc();
	} while (c == PEOA);
	return c;
}


static int
preadfd(void)
{
	int nr;
	char *buf =  parsefile->buf;
	parsefile->nextc = buf;

retry:
#ifndef SMALL
	if (parsefile->fd == 0 && el) {
		static const char *rl_cp;
		static int el_len;

		if (rl_cp == NULL) {
			struct stackmark smark;
			pushstackmark(&smark, stackblocksize());
			rl_cp = el_gets(el, &el_len);
			popstackmark(&smark);
		}
		if (rl_cp == NULL)
			nr = 0;
		else {
			nr = el_len;
			if (nr > IBUFSIZ - 1)
				nr = IBUFSIZ - 1;
			memcpy(buf, rl_cp, nr);
			if (nr != el_len) {
				el_len -= nr;
				rl_cp += nr;
			} else
				rl_cp = 0;
		}

	} else
#endif
		nr = read(parsefile->fd, buf, IBUFSIZ - 1);


	if (nr < 0) {
		if (errno == EINTR)
			goto retry;
		if (parsefile->fd == 0 && errno == EWOULDBLOCK) {
			int flags = fcntl(0, F_GETFL, 0);
			if (flags >= 0 && flags & O_NONBLOCK) {
				flags &=~ O_NONBLOCK;
				if (fcntl(0, F_SETFL, flags) >= 0) {
					out2str("sh: turning off NDELAY mode\n");
					goto retry;
				}
			}
		}
	}
	return nr;
}

/*
 * Refill the input buffer and return the next input character:
 *
 * 1) If a string was pushed back on the input, pop it;
 * 2) If an EOF was pushed back (parsenleft == EOF_NLEFT) or we are reading
 *    from a string so we can't refill the buffer, return EOF.
 * 3) If the is more stuff in this buffer, use it else call read to fill it.
 * 4) Process input up to the next newline, deleting nul characters.
 */

static int preadbuffer(void)
{
	char *q;
	int more;
#ifndef SMALL
	int something;
#endif
	char savec;

	if (unlikely(parsefile->strpush)) {
		if (
			parsefile->nleft == -1 &&
			parsefile->strpush->ap &&
			parsefile->nextc[-1] != ' ' &&
			parsefile->nextc[-1] != '\t'
		) {
			return PEOA;
		}
		popstring();
		return pgetc();
	}
	if (unlikely(parsefile->nleft == EOF_NLEFT ||
		     parsefile->buf == NULL))
		return PEOF;
	flushall();

	more = parsefile->lleft;
	if (more <= 0) {
again:
		if ((more = preadfd()) <= 0) {
			parsefile->lleft = parsefile->nleft = EOF_NLEFT;
			return PEOF;
		}
	}

	q = parsefile->nextc;

	/* delete nul characters */
#ifndef SMALL
	something = 0;
#endif
	for (;;) {
		int c;

		more--;
		c = *q;

		if (!c)
			memmove(q, q + 1, more);
		else {
			q++;

			if (c == '\n') {
				parsefile->nleft = q - parsefile->nextc - 1;
				break;
			}

#ifndef SMALL
			switch (c) {
			default:
				something = 1;
				/* fall through */
			case '\t':
			case ' ':
				break;
			}
#endif
		}

		if (more <= 0) {
			parsefile->nleft = q - parsefile->nextc - 1;
			if (parsefile->nleft < 0)
				goto again;
			break;
		}
	}
	parsefile->lleft = more;

	savec = *q;
	*q = '\0';

#ifndef SMALL
	if (parsefile->fd == 0 && hist && something) {
		HistEvent he;
		INTOFF;
		history(hist, &he, whichprompt == 1? H_ENTER : H_APPEND,
			parsefile->nextc);
		INTON;
	}
#endif

	if (vflag) {
		out2str(parsefile->nextc);
#ifdef FLUSHERR
		flushout(out2);
#endif
	}

	*q = savec;

	return (signed char)*parsefile->nextc++;
}

/*
 * Undo a call to pgetc.  Only two characters may be pushed back.
 * PEOF may be pushed back.
 */

void
pungetc(void)
{
	parsefile->unget++;
}

/*
 * Push a string back onto the input at this current parsefile level.
 * We handle aliases this way.
 */
void
pushstring(char *s, void *ap)
{
	struct strpush *sp;
	size_t len;

	len = strlen(s);
	INTOFF;
/*dprintf("*** calling pushstring: %s, %d\n", s, len);*/
	if (parsefile->strpush) {
		sp = ckmalloc(sizeof (struct strpush));
		sp->prev = parsefile->strpush;
		parsefile->strpush = sp;
	} else
		sp = parsefile->strpush = &(parsefile->basestrpush);
	sp->prevstring = parsefile->nextc;
	sp->prevnleft = parsefile->nleft;
	sp->unget = parsefile->unget;
	memcpy(sp->lastc, parsefile->lastc, sizeof(sp->lastc));
	sp->ap = (struct alias *)ap;
	if (ap) {
		((struct alias *)ap)->flag |= ALIASINUSE;
		sp->string = s;
	}
	parsefile->nextc = s;
	parsefile->nleft = len;
	parsefile->unget = 0;
	INTON;
}

void
popstring(void)
{
	struct strpush *sp = parsefile->strpush;

	INTOFF;
	if (sp->ap) {
		if (parsefile->nextc[-1] == ' ' ||
		    parsefile->nextc[-1] == '\t') {
			checkkwd |= CHKALIAS;
		}
		if (sp->string != sp->ap->val) {
			ckfree(sp->string);
		}
		sp->ap->flag &= ~ALIASINUSE;
		if (sp->ap->flag & ALIASDEAD) {
			unalias(sp->ap->name);
		}
	}
	parsefile->nextc = sp->prevstring;
	parsefile->nleft = sp->prevnleft;
	parsefile->unget = sp->unget;
	memcpy(parsefile->lastc, sp->lastc, sizeof(sp->lastc));
/*dprintf("*** calling popstring: restoring to '%s'\n", parsenextc);*/
	parsefile->strpush = sp->prev;
	if (sp != &(parsefile->basestrpush))
		ckfree(sp);
	INTON;
}

/*
 * Set the input to take input from a file.  If push is set, push the
 * old input onto the stack first.
 */

int
setinputfile(const char *fname, int flags)
{
	int fd;

	INTOFF;
	if ((fd = open(fname, O_RDONLY)) < 0) {
		if (flags & INPUT_NOFILE_OK)
			goto out;
		exitstatus = 127;
		exerror(EXERROR, "Can't open %s", fname);
	}
	if (fd < 10)
		fd = savefd(fd, fd);
	setinputfd(fd, flags & INPUT_PUSH_FILE);
out:
	INTON;
	return fd;
}


/*
 * Like setinputfile, but takes an open file descriptor.  Call this with
 * interrupts off.
 */

/* 
static void
*/
void // libdash
setinputfd(int fd, int push)
{
	if (push) {
		pushfile();
		parsefile->buf = 0;
	}
	parsefile->fd = fd;
	if (parsefile->buf == NULL)
		parsefile->buf = ckmalloc(IBUFSIZ);
	parsefile->lleft = parsefile->nleft = 0;
	plinno = 1;
}


/*
 * Like setinputfile, but takes input from a string.
 */

void
setinputstring(char *string)
{
	INTOFF;
	pushfile();
	parsefile->nextc = string;
	parsefile->nleft = strlen(string);
	parsefile->buf = NULL;
	plinno = 1;
	INTON;
}



/*
 * To handle the "." command, a stack of input files is used.  Pushfile
 * adds a new entry to the stack and popfile restores the previous level.
 */

STATIC void
pushfile(void)
{
	struct parsefile *pf;

	pf = (struct parsefile *)ckmalloc(sizeof (struct parsefile));
	pf->prev = parsefile;
	pf->fd = -1;
	pf->strpush = NULL;
	pf->basestrpush.prev = NULL;
	pf->unget = 0;
	parsefile = pf;
}


void
popfile(void)
{
	struct parsefile *pf = parsefile;

	INTOFF;
	if (pf->fd >= 0)
		close(pf->fd);
	if (pf->buf)
		ckfree(pf->buf);
	while (pf->strpush)
		popstring();
	parsefile = pf->prev;
	ckfree(pf);
	INTON;
}


void unwindfiles(struct parsefile *stop)
{
	while (parsefile != stop)
		popfile();
}


/*
 * Return to top level.
 */

void
popallfiles(void)
{
	unwindfiles(&basepf);
}



/*
 * Close the file(s) that the shell is reading commands from.  Called
 * after a fork is done.
 */

void
closescript(void)
{
	popallfiles();
	if (parsefile->fd > 0) {
		close(parsefile->fd);
		parsefile->fd = 0;
	}
}
