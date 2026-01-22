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

#include <sys/types.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <dirent.h>
#include <unistd.h>
#ifdef HAVE_GETPWNAM
#include <pwd.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include <limits.h>
#include <string.h>
#ifdef HAVE_FNMATCH
#include <fnmatch.h>
#endif
#ifdef HAVE_GLOB
#include <glob.h>
#endif
#include <ctype.h>
#include <stdbool.h>

/*
 * Routines to expand arguments to commands.  We have to deal with
 * backquotes, shell variables, and file metacharacters.
 */

#include "shell.h"
#include "main.h"
#include "nodes.h"
#include "eval.h"
#include "expand.h"
#include "syntax.h"
#include "parser.h"
#include "jobs.h"
#include "options.h"
#include "var.h"
#include "output.h"
#include "memalloc.h"
#include "error.h"
#include "mystring.h"
#include "show.h"
#include "system.h"

/*
 * _rmescape() flags
 */
#define RMESCAPE_ALLOC	0x1	/* Allocate a new string */
#define RMESCAPE_GLOB	0x2	/* Add backslashes for glob */
#define RMESCAPE_GROW	0x8	/* Grow strings instead of stalloc */
#define RMESCAPE_HEAP	0x10	/* Malloc strings instead of stalloc */

/* Add CTLESC when necessary. */
#define QUOTES_ESC	(EXP_FULL | EXP_CASE)

/*
 * Structure specifying which parts of the string should be searched
 * for IFS characters.
 */

struct ifsregion {
	struct ifsregion *next;	/* next region in list */
	int begoff;		/* offset of start of region */
	int endoff;		/* offset of end of region */
	int nulonly;		/* search for nul bytes only */
};

/* output of current string */
static char *expdest;
/* list of back quote expressions */
static struct nodelist *argbackq;
/* first struct in list of ifs regions */
static struct ifsregion ifsfirst;
/* last struct in list */
static struct ifsregion *ifslastp;
/* holds expanded arg list */
static struct arglist exparg;

static char *argstr(char *p, int flag);
static char *exptilde(char *startp, int flag);
static char *expari(char *start, int flag);
STATIC void expbackq(union node *, int);
STATIC char *evalvar(char *, int);
static size_t strtodest(const char *p, int flags);
static size_t memtodest(const char *p, size_t len, int flags);
STATIC ssize_t varvalue(char *, int, int, int);
STATIC void expandmeta(struct strlist *, int);
#ifdef HAVE_GLOB
STATIC void addglob(const glob_t *);
#else
STATIC void expmeta(char *, unsigned, unsigned);
STATIC struct strlist *expsort(struct strlist *);
STATIC struct strlist *msort(struct strlist *, int);
#endif
STATIC void addfname(char *);
STATIC int patmatch(char *, const char *);
#ifndef HAVE_FNMATCH
STATIC int pmatch(const char *, const char *);
#else
#define pmatch(a, b) !fnmatch((a), (b), 0)
#endif
static size_t cvtnum(intmax_t num, int flags);
STATIC size_t esclen(const char *, const char *);
STATIC char *scanleft(char *, char *, char *, char *, int, int);
STATIC char *scanright(char *, char *, char *, char *, int, int);
STATIC void varunset(const char *, const char *, const char *, int)
	__attribute__((__noreturn__));


/*
 * Prepare a pattern for a glob(3) call.
 *
 * Returns an stalloced string.
 */

STATIC inline char *
preglob(const char *pattern, int flag) {
	flag |= RMESCAPE_GLOB;
	return _rmescapes((char *)pattern, flag);
}


STATIC size_t
esclen(const char *start, const char *p) {
	size_t esc = 0;

	while (p > start && *--p == (char)CTLESC) {
		esc++;
	}
	return esc;
}


static inline const char *getpwhome(const char *name)
{
#ifdef HAVE_GETPWNAM
	struct passwd *pw = getpwnam(name);
	return pw ? pw->pw_dir : 0;
#else
	return 0;
#endif
}


/*
 * Perform variable substitution and command substitution on an argument,
 * placing the resulting list of arguments in arglist.  If EXP_FULL is true,
 * perform splitting and file name expansion.  When arglist is NULL, perform
 * here document expansion.
 */

void
expandarg(union node *arg, struct arglist *arglist, int flag)
{
	struct strlist *sp;
	char *p;

	argbackq = arg->narg.backquote;
	STARTSTACKSTR(expdest);
	argstr(arg->narg.text, flag);
	if (arglist == NULL) {
		/* here document expanded */
		goto out;
	}
	p = grabstackstr(expdest);
	exparg.lastp = &exparg.list;
	/*
	 * TODO - EXP_REDIR
	 */
	if (flag & EXP_FULL) {
		ifsbreakup(p, -1, &exparg);
		*exparg.lastp = NULL;
		exparg.lastp = &exparg.list;
		expandmeta(exparg.list, flag);
	} else {
		sp = (struct strlist *)stalloc(sizeof (struct strlist));
		sp->text = p;
		*exparg.lastp = sp;
		exparg.lastp = &sp->next;
	}
	*exparg.lastp = NULL;
	if (exparg.list) {
		*arglist->lastp = exparg.list;
		arglist->lastp = exparg.lastp;
	}

out:
	ifsfree();
}



/*
 * Perform variable and command substitution.  If EXP_FULL is set, output CTLESC
 * characters to allow for further processing.  Otherwise treat
 * $@ like $* since no splitting will be performed.
 */

static char *argstr(char *p, int flag)
{
	static const char spclchars[] = {
		'=',
		':',
		CTLQUOTEMARK,
		CTLENDVAR,
		CTLESC,
		CTLVAR,
		CTLBACKQ,
		CTLARI,
		CTLENDARI,
		0
	};
	const char *reject = spclchars;
	int c;
	int breakall = (flag & (EXP_WORD | EXP_QUOTED)) == EXP_WORD;
	int inquotes;
	size_t length;
	int startloc;

	reject += !!(flag & EXP_VARTILDE2);
	reject += flag & EXP_VARTILDE ? 0 : 2;
	inquotes = 0;
	length = 0;
	if (flag & EXP_TILDE) {
		flag &= ~EXP_TILDE;
tilde:
		if (*p == '~')
			p = exptilde(p, flag);
	}
start:
	startloc = expdest - (char *)stackblock();
	for (;;) {
		int end;

		length += strcspn(p + length, reject);
		end = 0;
		c = (signed char)p[length];
		if (!(c & 0x80) || c == CTLENDARI || c == CTLENDVAR) {
			/*
			 * c == '=' || c == ':' || c == '\0' ||
			 * c == CTLENDARI || c == CTLENDVAR
			 */
			length++;
			/* c == '\0' || c == CTLENDARI || c == CTLENDVAR */
			end = !!((c - 1) & 0x80);
		}
		if (length > 0 && !(flag & EXP_DISCARD)) {
			int newloc;
			char *q;

			q = stnputs(p, length, expdest);
			q[-1] &= end - 1;
			expdest = q - (flag & EXP_WORD ? end : 0);
			newloc = expdest - (char *)stackblock() - end;
			if (breakall && !inquotes && newloc > startloc) {
				recordregion(startloc, newloc, 0);
			}
			startloc = newloc;
		}
		p += length + 1;
		length = 0;

		if (end)
			break;

		switch (c) {
		case '=':
			flag |= EXP_VARTILDE2;
			reject++;
			/* fall through */
		case ':':
			/*
			 * sort of a hack - expand tildes in variable
			 * assignments (after the first '=' and after ':'s).
			 */
			if (*--p == '~') {
				goto tilde;
			}
			continue;
		case CTLQUOTEMARK:
			/* "$@" syntax adherence hack */
			if (!inquotes && !memcmp(p, dolatstr + 1,
						 DOLATSTRLEN - 1)) {
				p = evalvar(p + 1, flag | EXP_QUOTED) + 1;
				goto start;
			}
			inquotes ^= EXP_QUOTED;
addquote:
			if (flag & QUOTES_ESC) {
				p--;
				length++;
				startloc++;
			}
			break;
		case CTLESC:
			startloc++;
			length++;
			goto addquote;
		case CTLVAR:
			p = evalvar(p, flag | inquotes);
			goto start;
		case CTLBACKQ:
			expbackq(argbackq->n, flag | inquotes);
			goto start;
		case CTLARI:
			p = expari(p, flag | inquotes);
			goto start;
		}
	}
	return p - 1;
}

static char *exptilde(char *startp, int flag)
{
	signed char c;
	char *name;
	const char *home;
	char *p;

	p = startp;
	name = p + 1;

	while ((c = *++p) != '\0') {
		switch(c) {
		case CTLESC:
			return (startp);
		case CTLQUOTEMARK:
			return (startp);
		case ':':
			if (flag & EXP_VARTILDE)
				goto done;
			break;
		case '/':
		case CTLENDVAR:
			goto done;
		}
	}
done:
	if (flag & EXP_DISCARD)
		goto out;
	*p = '\0';
	if (*name == '\0') {
		home = lookupvar(homestr);
	} else {
		home = getpwhome(name);
	}
	*p = c;
	if (!home)
		goto lose;
	strtodest(home, flag | EXP_QUOTED);
out:
	return (p);
lose:
	return (startp);
}


void 
removerecordregions(int endoff)
{
	if (ifslastp == NULL)
		return;

	if (ifsfirst.endoff > endoff) {
		while (ifsfirst.next != NULL) {
			struct ifsregion *ifsp;
			INTOFF;
			ifsp = ifsfirst.next->next;
			ckfree(ifsfirst.next);
			ifsfirst.next = ifsp;
			INTON;
		}
		if (ifsfirst.begoff > endoff)
			ifslastp = NULL;
		else {
			ifslastp = &ifsfirst;
			ifsfirst.endoff = endoff;
		}
		return;
	}
	
	ifslastp = &ifsfirst;
	while (ifslastp->next && ifslastp->next->begoff < endoff)
		ifslastp=ifslastp->next;
	while (ifslastp->next != NULL) {
		struct ifsregion *ifsp;
		INTOFF;
		ifsp = ifslastp->next->next;
		ckfree(ifslastp->next);
		ifslastp->next = ifsp;
		INTON;
	}
	if (ifslastp->endoff > endoff)
		ifslastp->endoff = endoff;
}


/*
 * Expand arithmetic expression.  Backup to start of expression,
 * evaluate, place result in (backed up) result, adjust string position.
 */
static char *expari(char *start, int flag)
{
	struct stackmark sm;
	int begoff;
	int endoff;
	int len;
	intmax_t result;
	char *p;

	p = stackblock();
	begoff = expdest - p;
	p = argstr(start, flag & EXP_DISCARD);

	if (flag & EXP_DISCARD)
		goto out;

	start = stackblock();
	endoff = expdest - start;
	start += begoff;
	STADJUST(start - expdest, expdest);

	removerecordregions(begoff);

	if (likely(flag & QUOTES_ESC))
		rmescapes(start);

	pushstackmark(&sm, endoff);
	result = arith(start);
	popstackmark(&sm);

	len = cvtnum(result, flag);

	if (likely(!(flag & EXP_QUOTED)))
		recordregion(begoff, begoff + len, 0);

out:
	return p;
}


/*
 * Expand stuff in backwards quotes.
 */

STATIC void
expbackq(union node *cmd, int flag)
{
	struct backcmd in;
	int i;
	char buf[128];
	char *p;
	char *dest;
	int startloc;
	struct stackmark smark;

	if (flag & EXP_DISCARD)
		goto out;

	INTOFF;
	startloc = expdest - (char *)stackblock();
	pushstackmark(&smark, startloc);
	evalbackcmd(cmd, (struct backcmd *) &in);
	popstackmark(&smark);

	p = in.buf;
	i = in.nleft;
	if (i == 0)
		goto read;
	for (;;) {
		memtodest(p, i, flag);
read:
		if (in.fd < 0)
			break;
		do {
			i = read(in.fd, buf, sizeof buf);
		} while (i < 0 && errno == EINTR);
		TRACE(("expbackq: read returns %d\n", i));
		if (i <= 0)
			break;
		p = buf;
	}

	if (in.buf)
		ckfree(in.buf);
	if (in.fd >= 0) {
		close(in.fd);
		back_exitstatus = waitforjob(in.jp);
	}
	INTON;

	/* Eat all trailing newlines */
	dest = expdest;
	for (; dest > (char *)stackblock() && dest[-1] == '\n';)
		STUNPUTC(dest);
	expdest = dest;

	if (!(flag & EXP_QUOTED))
		recordregion(startloc, dest - (char *)stackblock(), 0);
	TRACE(("evalbackq: size=%d: \"%.*s\"\n",
		(dest - (char *)stackblock()) - startloc,
		(dest - (char *)stackblock()) - startloc,
		stackblock() + startloc));

out:
	argbackq = argbackq->next;
}


STATIC char *
scanleft(
	char *startp, char *rmesc, char *rmescend, char *str, int quotes,
	int zero
) {
	char *loc;
	char *loc2;
	char c;

	loc = startp;
	loc2 = rmesc;
	do {
		int match;
		const char *s = loc2;
		c = *loc2;
		if (zero) {
			*loc2 = '\0';
			s = rmesc;
		}
		match = pmatch(str, s);
		*loc2 = c;
		if (match)
			return loc;
		if (quotes && *loc == (char)CTLESC)
			loc++;
		loc++;
		loc2++;
	} while (c);
	return 0;
}


STATIC char *
scanright(
	char *startp, char *rmesc, char *rmescend, char *str, int quotes,
	int zero
) {
	int esc = 0;
	char *loc;
	char *loc2;

	for (loc = str - 1, loc2 = rmescend; loc >= startp; loc2--) {
		int match;
		char c = *loc2;
		const char *s = loc2;
		if (zero) {
			*loc2 = '\0';
			s = rmesc;
		}
		match = pmatch(str, s);
		*loc2 = c;
		if (match)
			return loc;
		loc--;
		if (quotes) {
			if (--esc < 0) {
				esc = esclen(startp, loc);
			}
			if (esc % 2) {
				esc--;
				loc--;
			}
		}
	}
	return 0;
}

static char *subevalvar(char *start, char *str, int strloc, int startloc,
			int varflags, int flag)
{
	int subtype = varflags & VSTYPE;
	int quotes = flag & QUOTES_ESC;
	char *startp;
	char *loc;
	long amount;
	char *rmesc, *rmescend;
	int zero;
	char *(*scan)(char *, char *, char *, char *, int , int);
	char *p;

	p = argstr(start, (flag & EXP_DISCARD) | EXP_TILDE |
			  (str ?  0 : EXP_CASE));
	if (flag & EXP_DISCARD)
		return p;

	startp = stackblock() + startloc;

	switch (subtype) {
	case VSASSIGN:
		setvar(str, startp, 0);

		loc = startp;
		goto out;

	case VSQUESTION:
		varunset(start, str, startp, varflags);
		/* NOTREACHED */
	}

	subtype -= VSTRIMRIGHT;
#ifdef DEBUG
	if (subtype < 0 || subtype > 3)
		abort();
#endif

	rmesc = startp;
	rmescend = stackblock() + strloc;
	if (quotes) {
		rmesc = _rmescapes(startp, RMESCAPE_ALLOC | RMESCAPE_GROW);
		if (rmesc != startp) {
			rmescend = expdest;
			startp = stackblock() + startloc;
		}
	}
	rmescend--;
	str = stackblock() + strloc;
	preglob(str, 0);

	/* zero = subtype == VSTRIMLEFT || subtype == VSTRIMLEFTMAX */
	zero = subtype >> 1;
	/* VSTRIMLEFT/VSTRIMRIGHTMAX -> scanleft */
	scan = (subtype & 1) ^ zero ? scanleft : scanright;

	loc = scan(startp, rmesc, rmescend, str, quotes, zero);
	if (loc) {
		if (zero) {
			memmove(startp, loc, str - loc);
			loc = startp + (str - loc) - 1;
		}
		*loc = '\0';
	} else
		loc = str - 1;

out:
	amount = loc - expdest;
	STADJUST(amount, expdest);

	/* Remove any recorded regions beyond start of variable */
	removerecordregions(startloc);

	return p;
}


/*
 * Expand a variable, and return a pointer to the next character in the
 * input string.
 */
STATIC char *
evalvar(char *p, int flag)
{
	int subtype;
	int varflags;
	char *var;
	int patloc;
	int startloc;
	ssize_t varlen;
	int discard;
	int quoted;

	varflags = *p++;
	subtype = varflags & VSTYPE;

	quoted = flag & EXP_QUOTED;
	var = p;
	startloc = expdest - (char *)stackblock();
	p = strchr(p, '=') + 1;

again:
	varlen = varvalue(var, varflags, flag, quoted);
	if (varflags & VSNUL)
		varlen--;

	discard = varlen < 0 ? EXP_DISCARD : 0;

	switch (subtype) {
	case VSPLUS:
		discard ^= EXP_DISCARD;
		/* fall through */

	case 0:
	case VSMINUS:
		p = argstr(p, flag | EXP_TILDE | EXP_WORD |
			      (discard ^ EXP_DISCARD));
		goto record;

	case VSASSIGN:
	case VSQUESTION:
		p = subevalvar(p, var, 0, startloc, varflags,
			       (flag & ~QUOTES_ESC) |
			       (discard ^ EXP_DISCARD));

		if ((flag | ~discard) & EXP_DISCARD)
			goto record;

		varflags &= ~VSNUL;
		subtype = VSNORMAL;
		goto again;
	}

	if ((discard & ~flag) && uflag)
		varunset(p, var, 0, 0);

	if (subtype == VSLENGTH) {
		p++;
		if (flag & EXP_DISCARD)
			return p;
		cvtnum(varlen > 0 ? varlen : 0, flag);
		goto really_record;
	}

	if (subtype == VSNORMAL)
		goto record;

#ifdef DEBUG
	switch (subtype) {
	case VSTRIMLEFT:
	case VSTRIMLEFTMAX:
	case VSTRIMRIGHT:
	case VSTRIMRIGHTMAX:
		break;
	default:
		abort();
	}
#endif

	flag |= discard;
	if (!(flag & EXP_DISCARD)) {
		/*
		 * Terminate the string and start recording the pattern
		 * right after it
		 */
		STPUTC('\0', expdest);
	}

	patloc = expdest - (char *)stackblock();
	p = subevalvar(p, NULL, patloc, startloc, varflags, flag);

record:
	if ((flag | discard) & EXP_DISCARD)
		return p;

really_record:
	if (quoted) {
		quoted = *var == '@' && shellparam.nparam;
		if (!quoted)
			return p;
	}
	recordregion(startloc, expdest - (char *)stackblock(), quoted);
	return p;
}


/*
 * Put a string on the stack.
 */

static size_t memtodest(const char *p, size_t len, int flags)
{
	const char *syntax = flags & EXP_QUOTED ? DQSYNTAX : BASESYNTAX;
	char *q;
	char *s;

	if (unlikely(!len))
		return 0;

	q = makestrspace(len * 2, expdest);
	s = q;

	do {
		int c = (signed char)*p++;
		if (c) {
			if ((flags & QUOTES_ESC) &&
			    ((syntax[c] == CCTL) ||
			     (flags & EXP_QUOTED && syntax[c] == CBACK)))
				USTPUTC(CTLESC, q);
		} else if (!(flags & EXP_KEEPNUL))
			continue;
		USTPUTC(c, q);
	} while (--len);

	expdest = q;
	return q - s;
}


static size_t strtodest(const char *p, int flags)
{
	size_t len = strlen(p);
	memtodest(p, len, flags);
	return len;
}



/*
 * Add the value of a specialized variable to the stack string.
 */

STATIC ssize_t
varvalue(char *name, int varflags, int flags, int quoted)
{
	int num;
	char *p;
	int i;
	int sep;
	char sepc;
	char **ap;
	int subtype = varflags & VSTYPE;
	int discard = (subtype == VSPLUS || subtype == VSLENGTH) |
		      (flags & EXP_DISCARD);
	ssize_t len = 0;
	char c;

	if (!subtype) {
		if (discard)
			return -1;

		sh_error("Bad substitution");
	}

	flags |= EXP_KEEPNUL;
	flags &= discard ? ~QUOTES_ESC : ~0;
	sep = (flags & EXP_FULL) << CHAR_BIT;

	switch (*name) {
	case '$':
		num = rootpid;
		goto numvar;
	case '?':
		num = exitstatus;
		goto numvar;
	case '#':
		num = shellparam.nparam;
		goto numvar;
	case '!':
		num = backgndpid;
		if (num == 0)
			return -1;
numvar:
		len = cvtnum(num, flags);
		break;
	case '-':
		p = makestrspace(NOPTS, expdest);
		for (i = NOPTS - 1; i >= 0; i--) {
			if (optlist[i] && optletters[i]) {
				USTPUTC(optletters[i], p);
				len++;
			}
		}
		expdest = p;
		break;
	case '@':
		if (quoted && sep)
			goto param;
		/* fall through */
	case '*':
		/* We will set c to 0 or ~0 depending on whether
		 * we're doing field splitting.  We won't do field
		 * splitting if either we're quoted or sep is zero.
		 *
		 * Instead of testing (quoted || !sep) the following
		 * trick optimises away any branches by using the
		 * fact that EXP_QUOTED (which is the only bit that
		 * can be set in quoted) is the same as EXP_FULL <<
		 * CHAR_BIT (which is the only bit that can be set
		 * in sep).
		 */
#if EXP_QUOTED >> CHAR_BIT != EXP_FULL
#error The following two lines expect EXP_QUOTED == EXP_FULL << CHAR_BIT
#endif
		c = !((quoted | ~sep) & EXP_QUOTED) - 1;
		sep &= ~quoted;
		sep |= ifsset() ? (unsigned char)(c & ifsval()[0]) : ' ';
param:
		sepc = sep;
		if (!(ap = shellparam.p))
			return -1;
		while ((p = *ap++)) {
			len += strtodest(p, flags);

			if (*ap && sep) {
				len++;
				memtodest(&sepc, 1, flags);
			}
		}
		break;
	case '0':
	case '1':
	case '2':
	case '3':
	case '4':
	case '5':
	case '6':
	case '7':
	case '8':
	case '9':
		num = atoi(name);
		if (num < 0 || num > shellparam.nparam)
			return -1;
		p = num ? shellparam.p[num - 1] : arg0;
		goto value;
	default:
		p = lookupvar(name);
value:
		if (!p)
			return -1;

		len = strtodest(p, flags);
		break;
	}

	if (discard)
		STADJUST(-len, expdest);

	return len;
}



/*
 * Record the fact that we have to scan this region of the
 * string for IFS characters.
 */

void
recordregion(int start, int end, int nulonly)
{
	struct ifsregion *ifsp;

	if (ifslastp == NULL) {
		ifsp = &ifsfirst;
	} else {
		INTOFF;
		ifsp = (struct ifsregion *)ckmalloc(sizeof (struct ifsregion));
		ifsp->next = NULL;
		ifslastp->next = ifsp;
		INTON;
	}
	ifslastp = ifsp;
	ifslastp->begoff = start;
	ifslastp->endoff = end;
	ifslastp->nulonly = nulonly;
}



/*
 * Break the argument string into pieces based upon IFS and add the
 * strings to the argument list.  The regions of the string to be
 * searched for IFS characters have been stored by recordregion.
 * If maxargs is non-negative, at most maxargs arguments will be created, by
 * joining together the last arguments.
 */
void
ifsbreakup(char *string, int maxargs, struct arglist *arglist)
{
	struct ifsregion *ifsp;
	struct strlist *sp;
	char *start;
	char *p;
	char *q;
	char *r = NULL;
	const char *ifs, *realifs;
	int ifsspc;
	int nulonly;


	start = string;
	if (ifslastp != NULL) {
		ifsspc = 0;
		nulonly = 0;
		realifs = ifsset() ? ifsval() : defifs;
		ifsp = &ifsfirst;
		do {
			int afternul;

			p = string + ifsp->begoff;
			afternul = nulonly;
			nulonly = ifsp->nulonly;
			ifs = nulonly ? nullstr : realifs;
			ifsspc = 0;
			while (p < string + ifsp->endoff) {
				int c;
				bool isifs;
				bool isdefifs;

				q = p;
				c = *p++;
				if (c == (char)CTLESC)
					c = *p++;

				isifs = strchr(ifs, c);
				isdefifs = false;
				if (isifs)
					isdefifs = strchr(defifs, c);

				/* If only reading one more argument:
				 * If we have exactly one field,
				 * read that field without its terminator.
				 * If we have more than one field,
				 * read all fields including their terminators,
				 * except for trailing IFS whitespace.
				 *
				 * This means that if we have only IFS
				 * characters left, and at most one
				 * of them is non-whitespace, we stop
				 * reading here.
				 * Otherwise, we read all the remaining
				 * characters except for trailing
				 * IFS whitespace.
				 *
				 * In any case, r indicates the start
				 * of the characters to remove, or NULL
				 * if no characters should be removed.
				 */
				if (!maxargs) {
					if (isdefifs) {
						if (!r)
							r = q;
						continue;
					}

					if (!(isifs && ifsspc))
						r = NULL;

					ifsspc = 0;
					continue;
				}

				if (ifsspc) {
					if (isifs)
						q = p;

					start = q;

					if (isdefifs)
						continue;

					isifs = false;
				}

				if (isifs) {
					if (!(afternul || nulonly))
						ifsspc = isdefifs;
					/* Ignore IFS whitespace at start */
					if (q == start && ifsspc) {
						start = p;
						ifsspc = 0;
						continue;
					}
					if (maxargs > 0 && !--maxargs) {
						r = q;
						continue;
					}
					*q = '\0';
					sp = (struct strlist *)stalloc(sizeof *sp);
					sp->text = start;
					*arglist->lastp = sp;
					arglist->lastp = &sp->next;
					start = p;
					continue;
				}

				ifsspc = 0;
			}
		} while ((ifsp = ifsp->next) != NULL);
		if (nulonly)
			goto add;
	}

	if (r)
		*r = '\0';

	if (!*start)
		return;

add:
	sp = (struct strlist *)stalloc(sizeof *sp);
	sp->text = start;
	*arglist->lastp = sp;
	arglist->lastp = &sp->next;
}

void ifsfree(void)
{
	struct ifsregion *p = ifsfirst.next;

	if (!p)
		goto out;

	INTOFF;
	do {
		struct ifsregion *ifsp;
		ifsp = p->next;
		ckfree(p);
		p = ifsp;
	} while (p);
	ifsfirst.next = NULL;
	INTON;

out:
	ifslastp = NULL;
}



/*
 * Expand shell metacharacters.  At this point, the only control characters
 * should be escapes.  The results are stored in the list exparg.
 */

#ifdef HAVE_GLOB
STATIC void
expandmeta(str, flag)
	struct strlist *str;
	int flag;
{
	/* TODO - EXP_REDIR */

	while (str) {
		const char *p;
		glob_t pglob;
		int i;

		if (fflag)
			goto nometa;
		INTOFF;
		p = preglob(str->text, RMESCAPE_ALLOC | RMESCAPE_HEAP);
		i = glob(p, GLOB_NOMAGIC, 0, &pglob);
		if (p != str->text)
			ckfree(p);
		switch (i) {
		case 0:
			if ((pglob.gl_flags & (GLOB_NOMAGIC | GLOB_NOCHECK)) ==
			    (GLOB_NOMAGIC | GLOB_NOCHECK))
				goto nometa2;
			addglob(&pglob);
			globfree(&pglob);
			INTON;
			break;
		case GLOB_NOMATCH:
nometa2:
			globfree(&pglob);
			INTON;
nometa:
			*exparg.lastp = str;
			rmescapes(str->text);
			exparg.lastp = &str->next;
			break;
		default:	/* GLOB_NOSPACE */
			sh_error("Out of space");
		}
		str = str->next;
	}
}


/*
 * Add the result of glob(3) to the list.
 */

STATIC void
addglob(pglob)
	const glob_t *pglob;
{
	char **p = pglob->gl_pathv;

	do {
		addfname(*p);
	} while (*++p);
}


#else	/* HAVE_GLOB */
STATIC char *expdir;
STATIC unsigned expdir_max;


STATIC void
expandmeta(struct strlist *str, int flag)
{
	static const char metachars[] = {
		'*', '?', '[', 0
	};
	/* TODO - EXP_REDIR */

	while (str) {
		struct strlist **savelastp;
		struct strlist *sp;
		char *p;
		unsigned len;

		if (fflag)
			goto nometa;
		if (!strpbrk(str->text, metachars))
			goto nometa;
		savelastp = exparg.lastp;

		INTOFF;
		p = preglob(str->text, RMESCAPE_ALLOC | RMESCAPE_HEAP);
		len = strlen(p);
		expdir_max = len + PATH_MAX;
		expdir = ckmalloc(expdir_max);

		expmeta(p, len, 0);
		ckfree(expdir);
		if (p != str->text)
			ckfree(p);
		INTON;
		if (exparg.lastp == savelastp) {
			/*
			 * no matches
			 */
nometa:
			*exparg.lastp = str;
			rmescapes(str->text);
			exparg.lastp = &str->next;
		} else {
			*exparg.lastp = NULL;
			*savelastp = sp = expsort(*savelastp);
			while (sp->next != NULL)
				sp = sp->next;
			exparg.lastp = &sp->next;
		}
		str = str->next;
	}
}


/*
 * Do metacharacter (i.e. *, ?, [...]) expansion.
 */

STATIC void
expmeta(char *name, unsigned name_len, unsigned expdir_len)
{
	char *enddir = expdir + expdir_len;
	char *p;
	const char *cp;
	char *start;
	char *endname;
	int metaflag;
	struct stat64 statb;
	DIR *dirp;
	struct dirent *dp;
	int atend;
	int matchdot;
	int esc;

	metaflag = 0;
	start = name;
	for (p = name; esc = 0, *p; p += esc + 1) {
		if (*p == '*' || *p == '?')
			metaflag = 1;
		else if (*p == '[') {
			char *q = p + 1;
			if (*q == '!')
				q++;
			for (;;) {
				if (*q == '\\')
					q++;
				if (*q == '/' || *q == '\0')
					break;
				if (*++q == ']') {
					metaflag = 1;
					break;
				}
			}
		} else {
			if (*p == '\\' && p[1])
				esc++;
			if (p[esc] == '/') {
				if (metaflag)
					break;
				start = p + esc + 1;
			}
		}
	}
	if (metaflag == 0) {	/* we've reached the end of the file name */
		if (!expdir_len)
			return;
		p = name;
		do {
			if (*p == '\\' && p[1])
				p++;
			*enddir++ = *p;
		} while (*p++);
		if (lstat64(expdir, &statb) >= 0)
			addfname(expdir);
		return;
	}
	endname = p;
	if (name < start) {
		p = name;
		do {
			if (*p == '\\' && p[1])
				p++;
			*enddir++ = *p++;
		} while (p < start);
	}
	*enddir = 0;
	cp = expdir;
	expdir_len = enddir - cp;
	if (!expdir_len)
		cp = ".";
	if ((dirp = opendir(cp)) == NULL)
		return;
	if (*endname == 0) {
		atend = 1;
	} else {
		atend = 0;
		*endname = '\0';
		endname += esc + 1;
	}
	name_len -= endname - name;
	matchdot = 0;
	p = start;
	if (*p == '\\')
		p++;
	if (*p == '.')
		matchdot++;
	while (! int_pending() && (dp = readdir(dirp)) != NULL) {
		if (dp->d_name[0] == '.' && ! matchdot)
			continue;
		if (pmatch(start, dp->d_name)) {
			if (atend) {
				scopy(dp->d_name, enddir);
				addfname(expdir);
			} else {
				unsigned offset;
				unsigned len;

				p = stpcpy(enddir, dp->d_name);
				*p = '/';

				offset = p - expdir + 1;
				len = offset + name_len + NAME_MAX;
				if (len > expdir_max) {
					len += PATH_MAX;
					expdir = ckrealloc(expdir, len);
					expdir_max = len;
				}

				expmeta(endname, name_len, offset);
				enddir = expdir + expdir_len;
			}
		}
	}
	closedir(dirp);
	if (! atend)
		endname[-esc - 1] = esc ? '\\' : '/';
}
#endif	/* HAVE_GLOB */


/*
 * Add a file name to the list.
 */

STATIC void
addfname(char *name)
{
	struct strlist *sp;

	sp = (struct strlist *)stalloc(sizeof *sp);
	sp->text = sstrdup(name);
	*exparg.lastp = sp;
	exparg.lastp = &sp->next;
}


#ifndef HAVE_GLOB
/*
 * Sort the results of file name expansion.  It calculates the number of
 * strings to sort and then calls msort (short for merge sort) to do the
 * work.
 */

STATIC struct strlist *
expsort(struct strlist *str)
{
	int len;
	struct strlist *sp;

	len = 0;
	for (sp = str ; sp ; sp = sp->next)
		len++;
	return msort(str, len);
}


STATIC struct strlist *
msort(struct strlist *list, int len)
{
	struct strlist *p, *q = NULL;
	struct strlist **lpp;
	int half;
	int n;

	if (len <= 1)
		return list;
	half = len >> 1;
	p = list;
	for (n = half ; --n >= 0 ; ) {
		q = p;
		p = p->next;
	}
	q->next = NULL;			/* terminate first half of list */
	q = msort(list, half);		/* sort first half of list */
	p = msort(p, len - half);		/* sort second half */
	lpp = &list;
	for (;;) {
		if (strcmp(p->text, q->text) < 0) {
			*lpp = p;
			lpp = &p->next;
			if ((p = *lpp) == NULL) {
				*lpp = q;
				break;
			}
		} else {
			*lpp = q;
			lpp = &q->next;
			if ((q = *lpp) == NULL) {
				*lpp = p;
				break;
			}
		}
	}
	return list;
}
#endif


/*
 * Returns true if the pattern matches the string.
 */

STATIC inline int
patmatch(char *pattern, const char *string)
{
	return pmatch(preglob(pattern, 0), string);
}


#ifndef HAVE_FNMATCH
STATIC int ccmatch(const char *p, int chr, const char **r)
{
	static const struct class {
		char name[10];
		int (*fn)(int);
	} classes[] = {
		{ .name = ":alnum:]", .fn = isalnum },
		{ .name = ":cntrl:]", .fn = iscntrl },
		{ .name = ":lower:]", .fn = islower },
		{ .name = ":space:]", .fn = isspace },
		{ .name = ":alpha:]", .fn = isalpha },
		{ .name = ":digit:]", .fn = isdigit },
		{ .name = ":print:]", .fn = isprint },
		{ .name = ":upper:]", .fn = isupper },
		{ .name = ":blank:]", .fn = isblank },
		{ .name = ":graph:]", .fn = isgraph },
		{ .name = ":punct:]", .fn = ispunct },
		{ .name = ":xdigit:]", .fn = isxdigit },
	};
	const struct class *class, *end;

	end = classes + sizeof(classes) / sizeof(classes[0]);
	for (class = classes; class < end; class++) {
		const char *q;

		q = prefix(p, class->name);
		if (!q)
			continue;
		*r = q;
		return class->fn(chr);
	}

	*r = 0;
	return 0;
}

STATIC int
pmatch(const char *pattern, const char *string)
{
	const char *p, *q;
	char c;

	p = pattern;
	q = string;
	for (;;) {
		switch (c = *p++) {
		case '\0':
			goto breakloop;
		case '\\':
			if (*p) {
				c = *p++;
			}
			goto dft;
		case '?':
			if (*q++ == '\0')
				return 0;
			break;
		case '*':
			c = *p;
			while (c == '*')
				c = *++p;
			if (c != '\\' && c != '?' && c != '*' && c != '[') {
				while (*q != c) {
					if (*q == '\0')
						return 0;
					q++;
				}
			}
			do {
				if (pmatch(p, q))
					return 1;
			} while (*q++ != '\0');
			return 0;
		case '[': {
			const char *startp;
			int invert, found;
			char chr;

			startp = p;
			invert = 0;
			if (*p == '!') {
				invert++;
				p++;
			}
			found = 0;
			chr = *q;
			if (chr == '\0')
				return 0;
			c = *p++;
			do {
				if (!c) {
					p = startp;
					c = '[';
					goto dft;
				}
				if (c == '[') {
					const char *r;

					found |= !!ccmatch(p, chr, &r);
					if (r) {
						p = r;
						continue;
					}
				} else if (c == '\\')
					c = *p++;
				if (*p == '-' && p[1] != ']') {
					p++;
					if (*p == '\\')
						p++;
					if (chr >= c && chr <= *p)
						found = 1;
					p++;
				} else {
					if (chr == c)
						found = 1;
				}
			} while ((c = *p++) != ']');
			if (found == invert)
				return 0;
			q++;
			break;
		}
dft:	        default:
			if (*q++ != c)
				return 0;
			break;
		}
	}
breakloop:
	if (*q != '\0')
		return 0;
	return 1;
}
#endif



/*
 * Remove any CTLESC characters from a string.
 */

char *
_rmescapes(char *str, int flag)
{
	char *p, *q, *r;
	int notescaped;
	int globbing;

	p = strpbrk(str, qchars);
	if (!p) {
		return str;
	}
	q = p;
	r = str;
	if (flag & RMESCAPE_ALLOC) {
		size_t len = p - str;
		size_t fulllen = len + strlen(p) + 1;

		if (flag & RMESCAPE_GROW) {
			int strloc = str - (char *)stackblock();

			r = makestrspace(fulllen, expdest);
			str = (char *)stackblock() + strloc;
			p = str + len;
		} else if (flag & RMESCAPE_HEAP) {
			r = ckmalloc(fulllen);
		} else {
			r = stalloc(fulllen);
		}
		q = r;
		if (len > 0) {
			q = mempcpy(q, str, len);
		}
	}
	globbing = flag & RMESCAPE_GLOB;
	notescaped = globbing;
	while (*p) {
		if (*p == (char)CTLQUOTEMARK) {
			p++;
			notescaped = globbing;
			continue;
		}
		if (*p == '\\') {
			/* naked back slash */
			notescaped = 0;
			goto copy;
		}
		if (*p == (char)CTLESC) {
			p++;
			if (notescaped)
				*q++ = '\\';
		}
		notescaped = globbing;
copy:
		*q++ = *p++;
	}
	*q = '\0';
	if (flag & RMESCAPE_GROW) {
		expdest = r;
		STADJUST(q - r + 1, expdest);
	}
	return r;
}



/*
 * See if a pattern matches in a case statement.
 */

int
casematch(union node *pattern, char *val)
{
	struct stackmark smark;
	int result;

	setstackmark(&smark);
	argbackq = pattern->narg.backquote;
	STARTSTACKSTR(expdest);
	argstr(pattern->narg.text, EXP_TILDE | EXP_CASE);
	ifsfree();
	result = patmatch(stackblock(), val);
	popstackmark(&smark);
	return result;
}

/*
 * Our own itoa().
 */

static size_t cvtnum(intmax_t num, int flags)
{
	int len = max_int_length(sizeof(num));
	char buf[len];

	len = fmtstr(buf, len, "%" PRIdMAX, num);
	return memtodest(buf, len, flags);
}

STATIC void
varunset(const char *end, const char *var, const char *umsg, int varflags)
{
	const char *msg;
	const char *tail;

	tail = nullstr;
	msg = "parameter not set";
	if (umsg) {
		if (*end == (char)CTLENDVAR) {
			if (varflags & VSNUL)
				tail = " or null";
		} else
			msg = umsg;
	}
	sh_error("%.*s: %s%s", end - var - 1, var, msg, tail);
}

#ifdef mkinit

INCLUDE "expand.h"

EXITRESET {
	ifsfree();
}

#endif
