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

#if HAVE_ALLOCA_H
#include <alloca.h>
#endif

#include <stdlib.h>

#include "shell.h"
#include "parser.h"
#include "nodes.h"
#include "expand.h"	/* defines rmescapes() */
#include "exec.h"	/* defines find_builtin() */
#include "syntax.h"
#include "options.h"
#include "input.h"
#include "output.h"
#include "var.h"
#include "error.h"
#include "memalloc.h"
#include "init.h"       /* defines reset() */ // libdash
#include "mystring.h"
#include "alias.h"
#include "show.h"
#include "builtins.h"
#include "system.h"
#ifndef SMALL
#include "myhistedit.h"
#endif

/*
 * Shell command parser.
 */

/* values returned by readtoken */
#include "token_vars.h"



/* Used by expandstr to get here-doc like behaviour. */
#define FAKEEOFMARK (char *)1



struct heredoc {
	struct heredoc *next;	/* next here document in list */
	union node *here;		/* redirection node */
	char *eofmark;		/* string indicating end of input */
	int striptabs;		/* if set, strip leading tabs */
};

struct synstack {
	const char *syntax;
	struct synstack *prev;
	struct synstack *next;
	int innerdq;
	int varpushed;
	int dblquote;
	int varnest;		/* levels of variables expansion */
	int parenlevel;		/* levels of parens in arithmetic */
	int dqvarnest;		/* levels of variables expansion within double quotes */
};



struct heredoc *heredoclist;	/* list of here documents to read */
int doprompt;			/* if set, prompt the user */
int needprompt;			/* true if interactive and at start of line */
int lasttoken;			/* last token read */
int tokpushback;		/* last token pushed back */
char *wordtext;			/* text of last word returned by readtoken */
int checkkwd;
struct nodelist *backquotelist;
union node *redirnode;
struct heredoc *heredoc;
int quoteflag;			/* set if (part of) last token was quoted */


STATIC union node *list(int);
STATIC union node *andor(void);
STATIC union node *pipeline(void);
STATIC union node *command(void);
STATIC union node *simplecmd(void);
STATIC union node *makename(void);
STATIC void parsefname(void);
STATIC void parseheredoc(void);
STATIC int peektoken(void);
STATIC int readtoken(void);
STATIC int xxreadtoken(void);
STATIC int pgetc_eatbnl();
STATIC int readtoken1(int, char const *, char *, int);
STATIC void synexpect(int) __attribute__((__noreturn__));
STATIC void synerror(const char *) __attribute__((__noreturn__));
STATIC void setprompt(int);


int isassignment(const char *p)
{
	const char *q = endofname(p);
	if (p == q)
		return 0;
	return *q == '=';
}

static inline int realeofmark(const char *eofmark)
{
	return eofmark && eofmark != FAKEEOFMARK;
}


/*
 * Read and parse a command.  Returns NEOF on end of file.  (NULL is a
 * valid parse tree indicating a blank line.)
 */

union node *
parsecmd(int interact)
{
	tokpushback = 0;
	checkkwd = 0;
	heredoclist = 0;
	doprompt = interact;
	if (doprompt)
		setprompt(doprompt);
	needprompt = 0;

	return list(1);
}

// libdash
/* 2018-09-25 manually install a handler here so we can return an appropriate error code */
union node *
parsecmd_safe(int interact)
{
	struct jmploc jmploc;

	tokpushback = 0;
	checkkwd = 0;
	heredoclist = 0;
	doprompt = interact;
	if (doprompt)
		setprompt(doprompt);
	needprompt = 0;

        if (unlikely(setjmp(jmploc.loc))) {
          return NERR;
        }
        handler = &jmploc;

	return list(1);
}

STATIC union node *
list(int nlflag)
{
	union node *n1, *n2, *n3;
	int tok;

	n1 = NULL;
	for (;;) {
		switch (readtoken()) {
		case TNL:
			if (!(nlflag & 1))
				break;
			parseheredoc();
			return n1;

		case TEOF:
			if (!n1 && (nlflag & 1))
				n1 = NEOF;
			parseheredoc();
			tokpushback++;
			lasttoken = TEOF;
			return n1;
		}

		tokpushback++;
		checkkwd = CHKNL | CHKKWD | CHKALIAS;
		if (nlflag == 2 && tokendlist[peektoken()])
			return n1;
		nlflag |= 2;

		n2 = andor();
		tok = readtoken();
		if (tok == TBACKGND) {
			if (n2->type == NPIPE) {
				n2->npipe.backgnd = 1;
			} else {
				if (n2->type != NREDIR) {
					n3 = stalloc(sizeof(struct nredir));
					n3->nredir.n = n2;
					n3->nredir.redirect = NULL;
					n2 = n3;
				}
				n2->type = NBACKGND;
			}
		}
		if (n1 == NULL) {
			n1 = n2;
		}
		else {
			n3 = (union node *)stalloc(sizeof (struct nbinary));
			n3->type = NSEMI;
			n3->nbinary.ch1 = n1;
			n3->nbinary.ch2 = n2;
			n1 = n3;
		}
		switch (tok) {
		case TNL:
		case TEOF:
			tokpushback++;
			/* fall through */
		case TBACKGND:
		case TSEMI:
			break;
		default:
			if ((nlflag & 1))
				synexpect(-1);
			tokpushback++;
			return n1;
		}
	}
}



STATIC union node *
andor(void)
{
	union node *n1, *n2, *n3;
	int t;

	n1 = pipeline();
	for (;;) {
		if ((t = readtoken()) == TAND) {
			t = NAND;
		} else if (t == TOR) {
			t = NOR;
		} else {
			tokpushback++;
			return n1;
		}
		checkkwd = CHKNL | CHKKWD | CHKALIAS;
		n2 = pipeline();
		n3 = (union node *)stalloc(sizeof (struct nbinary));
		n3->type = t;
		n3->nbinary.ch1 = n1;
		n3->nbinary.ch2 = n2;
		n1 = n3;
	}
}



STATIC union node *
pipeline(void)
{
	union node *n1, *n2, *pipenode;
	struct nodelist *lp, *prev;
	int negate;

	negate = 0;
	TRACE(("pipeline: entered\n"));
	if (readtoken() == TNOT) {
		negate = !negate;
		checkkwd = CHKKWD | CHKALIAS;
	} else
		tokpushback++;
	n1 = command();
	if (readtoken() == TPIPE) {
		pipenode = (union node *)stalloc(sizeof (struct npipe));
		pipenode->type = NPIPE;
		pipenode->npipe.backgnd = 0;
		lp = (struct nodelist *)stalloc(sizeof (struct nodelist));
		pipenode->npipe.cmdlist = lp;
		lp->n = n1;
		do {
			prev = lp;
			lp = (struct nodelist *)stalloc(sizeof (struct nodelist));
			checkkwd = CHKNL | CHKKWD | CHKALIAS;
			lp->n = command();
			prev->next = lp;
		} while (readtoken() == TPIPE);
		lp->next = NULL;
		n1 = pipenode;
	}
	tokpushback++;
	if (negate) {
		n2 = (union node *)stalloc(sizeof (struct nnot));
		n2->type = NNOT;
		n2->nnot.com = n1;
		return n2;
	} else
		return n1;
}



STATIC union node *
command(void)
{
	union node *n1, *n2;
	union node *ap, **app;
	union node *cp, **cpp;
	union node *redir, **rpp;
	union node **rpp2;
	int t;
	int savelinno;

	redir = NULL;
	rpp2 = &redir;

	savelinno = plinno;

	switch (readtoken()) {
	default:
		synexpect(-1);
		/* NOTREACHED */
	case TIF:
		n1 = (union node *)stalloc(sizeof (struct nif));
		n1->type = NIF;
		n1->nif.test = list(0);
		if (readtoken() != TTHEN)
			synexpect(TTHEN);
		n1->nif.ifpart = list(0);
		n2 = n1;
		while (readtoken() == TELIF) {
			n2->nif.elsepart = (union node *)stalloc(sizeof (struct nif));
			n2 = n2->nif.elsepart;
			n2->type = NIF;
			n2->nif.test = list(0);
			if (readtoken() != TTHEN)
				synexpect(TTHEN);
			n2->nif.ifpart = list(0);
		}
		if (lasttoken == TELSE)
			n2->nif.elsepart = list(0);
		else {
			n2->nif.elsepart = NULL;
			tokpushback++;
		}
		t = TFI;
		break;
	case TWHILE:
	case TUNTIL: {
		int got;
		n1 = (union node *)stalloc(sizeof (struct nbinary));
		n1->type = (lasttoken == TWHILE)? NWHILE : NUNTIL;
		n1->nbinary.ch1 = list(0);
		if ((got=readtoken()) != TDO) {
TRACE(("expecting DO got %s %s\n", tokname[got], got == TWORD ? wordtext : ""));
			synexpect(TDO);
		}
		n1->nbinary.ch2 = list(0);
		t = TDONE;
		break;
	}
	case TFOR:
		if (readtoken() != TWORD || quoteflag || ! goodname(wordtext))
			synerror("Bad for loop variable");
		n1 = (union node *)stalloc(sizeof (struct nfor));
		n1->type = NFOR;
		n1->nfor.linno = savelinno;
		n1->nfor.var = wordtext;
		checkkwd = CHKNL | CHKKWD | CHKALIAS;
		if (readtoken() == TIN) {
			app = &ap;
			while (readtoken() == TWORD) {
				n2 = (union node *)stalloc(sizeof (struct narg));
				n2->type = NARG;
				n2->narg.text = wordtext;
				n2->narg.backquote = backquotelist;
				*app = n2;
				app = &n2->narg.next;
			}
			*app = NULL;
			n1->nfor.args = ap;
			if (lasttoken != TNL && lasttoken != TSEMI)
				synexpect(-1);
		} else {
			n2 = (union node *)stalloc(sizeof (struct narg));
			n2->type = NARG;
			n2->narg.text = (char *)dolatstr;
			n2->narg.backquote = NULL;
			n2->narg.next = NULL;
			n1->nfor.args = n2;
			/*
			 * Newline or semicolon here is optional (but note
			 * that the original Bourne shell only allowed NL).
			 */
			if (lasttoken != TSEMI)
				tokpushback++;
		}
		checkkwd = CHKNL | CHKKWD | CHKALIAS;
		if (readtoken() != TDO)
			synexpect(TDO);
		n1->nfor.body = list(0);
		t = TDONE;
		break;
	case TCASE:
		n1 = (union node *)stalloc(sizeof (struct ncase));
		n1->type = NCASE;
		n1->ncase.linno = savelinno;
		if (readtoken() != TWORD)
			synexpect(TWORD);
		n1->ncase.expr = n2 = (union node *)stalloc(sizeof (struct narg));
		n2->type = NARG;
		n2->narg.text = wordtext;
		n2->narg.backquote = backquotelist;
		n2->narg.next = NULL;
		checkkwd = CHKNL | CHKKWD | CHKALIAS;
		if (readtoken() != TIN)
			synexpect(TIN);
		cpp = &n1->ncase.cases;
next_case:
		checkkwd = CHKNL | CHKKWD;
		t = readtoken();
		while(t != TESAC) {
			if (lasttoken == TLP)
				readtoken();
			*cpp = cp = (union node *)stalloc(sizeof (struct nclist));
			cp->type = NCLIST;
			app = &cp->nclist.pattern;
			for (;;) {
				*app = ap = (union node *)stalloc(sizeof (struct narg));
				ap->type = NARG;
				ap->narg.text = wordtext;
				ap->narg.backquote = backquotelist;
				if (readtoken() != TPIPE)
					break;
				app = &ap->narg.next;
				readtoken();
			}
			ap->narg.next = NULL;
			if (lasttoken != TRP)
				synexpect(TRP);
			cp->nclist.body = list(2);

			cpp = &cp->nclist.next;

			checkkwd = CHKNL | CHKKWD;
			if ((t = readtoken()) != TESAC) {
				if (t != TENDCASE)
					synexpect(TENDCASE);
				else
					goto next_case;
			}
		}
		*cpp = NULL;
		goto redir;
	case TLP:
		n1 = (union node *)stalloc(sizeof (struct nredir));
		n1->type = NSUBSHELL;
		n1->nredir.linno = savelinno;
		n1->nredir.n = list(0);
		n1->nredir.redirect = NULL;
		t = TRP;
		break;
	case TBEGIN:
		n1 = list(0);
		t = TEND;
		break;
	case TWORD:
	case TREDIR:
// libdash
/* 2019-04-25 to allow for proper handling of empty aliases */
        case TNL:
		tokpushback++;
		return simplecmd();
	}

	if (readtoken() != t)
		synexpect(t);

redir:
	/* Now check for redirection which may follow command */
	checkkwd = CHKKWD | CHKALIAS;
	rpp = rpp2;
	while (readtoken() == TREDIR) {
		*rpp = n2 = redirnode;
		rpp = &n2->nfile.next;
		parsefname();
	}
	tokpushback++;
	*rpp = NULL;
	if (redir) {
		if (n1->type != NSUBSHELL) {
			n2 = (union node *)stalloc(sizeof (struct nredir));
			n2->type = NREDIR;
			n2->nredir.linno = savelinno;
			n2->nredir.n = n1;
			n1 = n2;
		}
		n1->nredir.redirect = redir;
	}

	return n1;
}


STATIC union node *
simplecmd(void) {
	union node *args, **app;
	union node *n = NULL;
	union node *vars, **vpp;
	union node **rpp, *redir;
	int savecheckkwd;
	int savelinno;

	args = NULL;
	app = &args;
	vars = NULL;
	vpp = &vars;
	redir = NULL;
	rpp = &redir;

	savecheckkwd = CHKALIAS;
	savelinno = plinno;
	for (;;) {
		checkkwd = savecheckkwd;
		switch (readtoken()) {
		case TWORD:
			n = (union node *)stalloc(sizeof (struct narg));
			n->type = NARG;
			n->narg.text = wordtext;
			n->narg.backquote = backquotelist;
			if (savecheckkwd && isassignment(wordtext)) {
				*vpp = n;
				vpp = &n->narg.next;
			} else {
				*app = n;
				app = &n->narg.next;
				savecheckkwd = 0;
			}
			break;
		case TREDIR:
			*rpp = n = redirnode;
			rpp = &n->nfile.next;
			parsefname();	/* read name of redirection file */
			break;
		case TLP:
			if (
				args && app == &args->narg.next &&
				!vars && !redir
			) {
				struct builtincmd *bcmd;
				const char *name;

				/* We have a function */
				if (readtoken() != TRP)
					synexpect(TRP);
				name = n->narg.text;
				if (
					!goodname(name) || (
						(bcmd = find_builtin(name)) &&
						bcmd->flags & BUILTIN_SPECIAL
					)
				)
					synerror("Bad function name");
				n->type = NDEFUN;
				checkkwd = CHKNL | CHKKWD | CHKALIAS;
				n->ndefun.text = n->narg.text;
				n->ndefun.linno = plinno;
				n->ndefun.body = command();
				return n;
			}
			/* fall through */
		default:
			tokpushback++;
			goto out;
		}
	}
out:
	*app = NULL;
	*vpp = NULL;
	*rpp = NULL;
	n = (union node *)stalloc(sizeof (struct ncmd));
	n->type = NCMD;
	n->ncmd.linno = savelinno;
	n->ncmd.args = args;
	n->ncmd.assign = vars;
	n->ncmd.redirect = redir;
	return n;
}

STATIC union node *
makename(void)
{
	union node *n;

	n = (union node *)stalloc(sizeof (struct narg));
	n->type = NARG;
	n->narg.next = NULL;
	n->narg.text = wordtext;
	n->narg.backquote = backquotelist;
	return n;
}

void fixredir(union node *n, const char *text, int err)
	{
	TRACE(("Fix redir %s %d\n", text, err));
	if (!err)
		n->ndup.vname = NULL;

	if (is_digit(text[0]) && text[1] == '\0')
		n->ndup.dupfd = digit_val(text[0]);
	else if (text[0] == '-' && text[1] == '\0')
		n->ndup.dupfd = -1;
	else {

		if (err)
			synerror("Bad fd number");
		else
			n->ndup.vname = makename();
	}
}


STATIC void
parsefname(void)
{
	union node *n = redirnode;

	if (n->type == NHERE)
		checkkwd = CHKEOFMARK;
	if (readtoken() != TWORD)
		synexpect(-1);
	if (n->type == NHERE) {
		struct heredoc *here = heredoc;
		struct heredoc *p;

		if (quoteflag == 0)
			n->type = NXHERE;
		TRACE(("Here document %d\n", n->type));
		rmescapes(wordtext);
		here->eofmark = wordtext;
		here->next = NULL;
		if (heredoclist == NULL)
			heredoclist = here;
		else {
			for (p = heredoclist ; p->next ; p = p->next);
			p->next = here;
		}
	} else if (n->type == NTOFD || n->type == NFROMFD) {
		fixredir(n, wordtext, 0);
	} else {
		n->nfile.fname = makename();
	}
}


/*
 * Input any here documents.
 */

STATIC void
parseheredoc(void)
{
	struct heredoc *here;
	union node *n;

	here = heredoclist;
	heredoclist = 0;

	while (here) {
		if (needprompt) {
			setprompt(2);
		}
		if (here->here->type == NHERE)
			readtoken1(pgetc(), SQSYNTAX, here->eofmark, here->striptabs);
		else
			readtoken1(pgetc_eatbnl(), DQSYNTAX, here->eofmark, here->striptabs);
		n = (union node *)stalloc(sizeof (struct narg));
		n->narg.type = NARG;
		n->narg.next = NULL;
		n->narg.text = wordtext;
		n->narg.backquote = backquotelist;
		here->here->nhere.doc = n;
		here = here->next;
	}
}

STATIC int
peektoken(void)
{
	int t;

	t = readtoken();
	tokpushback++;
	return (t);
}

STATIC int
readtoken(void)
{
	int t;
	int kwd = checkkwd;
#ifdef DEBUG
	int alreadyseen = tokpushback;
#endif

top:
	t = xxreadtoken();

	/*
	 * eat newlines
	 */
	if (kwd & CHKNL) {
		while (t == TNL) {
			parseheredoc();
			t = xxreadtoken();
		}
	}

// libdash
/* 2019-04-25 to handle empty aliases */
ignorenl:
	if (t != TWORD || quoteflag) {
		goto out;
	}

	/*
	 * check for keywords
	 */
	if (kwd & CHKKWD) {
		const char *const *pp;

		if ((pp = findkwd(wordtext))) {
			lasttoken = t = pp - parsekwd + KWDOFFSET;
			TRACE(("keyword %s recognized\n", tokname[t]));
			goto out;
		}
	}

	if (checkkwd & CHKALIAS) {
		struct alias *ap;
		if ((ap = lookupalias(wordtext, 1)) != NULL) {
// libdash
/* 2019-04-25 to handle empty aliases */
			if (*ap->val) {
				pushstring(ap->val, ap);
				goto top;
			} else {
				t = xxreadtoken();
				goto ignorenl;
   		        }
		}
	}
out:
	checkkwd = 0;
#ifdef DEBUG
	if (!alreadyseen)
	    TRACE(("token %s %s\n", tokname[t], t == TWORD ? wordtext : ""));
	else
	    TRACE(("reread token %s %s\n", tokname[t], t == TWORD ? wordtext : ""));
#endif
	return (t);
}

static void nlprompt(void)
{
	plinno++;
	if (doprompt)
		setprompt(2);
}

static void nlnoprompt(void)
{
	plinno++;
	needprompt = doprompt;
}


/*
 * Read the next input token.
 * If the token is a word, we set backquotelist to the list of cmds in
 *	backquotes.  We set quoteflag to true if any part of the word was
 *	quoted.
 * If the token is TREDIR, then we set redirnode to a structure containing
 *	the redirection.
 *
 * [Change comment:  here documents and internal procedures]
 * [Readtoken shouldn't have any arguments.  Perhaps we should make the
 *  word parsing code into a separate routine.  In this case, readtoken
 *  doesn't need to have any internal procedures, but parseword does.
 *  We could also make parseoperator in essence the main routine, and
 *  have parseword (readtoken1?) handle both words and redirection.]
 */

#define RETURN(token)	return lasttoken = token

STATIC int
xxreadtoken(void)
{
	int c;

	if (tokpushback) {
		tokpushback = 0;
		return lasttoken;
	}
	if (needprompt) {
		setprompt(2);
	}
	for (;;) {	/* until token or start of word found */
		c = pgetc_eatbnl();
		switch (c) {
		case ' ': case '\t':
		case PEOA:
			continue;
		case '#':
			while ((c = pgetc()) != '\n' && c != PEOF);
			pungetc();
			continue;
		case '\n':
			nlnoprompt();
			RETURN(TNL);
		case PEOF:
			RETURN(TEOF);
		case '&':
			if (pgetc_eatbnl() == '&')
				RETURN(TAND);
			pungetc();
			RETURN(TBACKGND);
		case '|':
			if (pgetc_eatbnl() == '|')
				RETURN(TOR);
			pungetc();
			RETURN(TPIPE);
		case ';':
			if (pgetc_eatbnl() == ';')
				RETURN(TENDCASE);
			pungetc();
			RETURN(TSEMI);
		case '(':
			RETURN(TLP);
		case ')':
			RETURN(TRP);
		}
		break;
	}
	return readtoken1(c, BASESYNTAX, (char *)NULL, 0);
#undef RETURN
}

static int pgetc_eatbnl(void)
{
	int c;

	while ((c = pgetc()) == '\\') {
		if (pgetc2() != '\n') {
			pungetc();
			break;
		}

		nlprompt();
	}

	return c;
}

static int pgetc_top(struct synstack *stack)
{
	return stack->syntax == SQSYNTAX ? pgetc() : pgetc_eatbnl();
}

static void synstack_push(struct synstack **stack, struct synstack *next,
			  const char *syntax)
{
	memset(next, 0, sizeof(*next));
	next->syntax = syntax;
	next->next = *stack;
	(*stack)->prev = next;
	*stack = next;
}

static void synstack_pop(struct synstack **stack)
{
	*stack = (*stack)->next;
}



/*
 * If eofmark is NULL, read a word or a redirection symbol.  If eofmark
 * is not NULL, read a here document.  In the latter case, eofmark is the
 * word which marks the end of the document and striptabs is true if
 * leading tabs should be stripped from the document.  The argument firstc
 * is the first character of the input token or document.
 *
 * Because C does not have internal subroutines, I have simulated them
 * using goto's to implement the subroutine linkage.  The following macros
 * will run code that appears at the end of readtoken1.
 */

#define CHECKEND()	{goto checkend; checkend_return:;}
#define PARSEREDIR()	{goto parseredir; parseredir_return:;}
#define PARSESUB()	{goto parsesub; parsesub_return:;}
#define PARSEBACKQOLD()	{oldstyle = 1; goto parsebackq; parsebackq_oldreturn:;}
#define PARSEBACKQNEW()	{oldstyle = 0; goto parsebackq; parsebackq_newreturn:;}
#define	PARSEARITH()	{goto parsearith; parsearith_return:;}

STATIC int
readtoken1(int firstc, char const *syntax, char *eofmark, int striptabs)
{
	int c = firstc;
	char *out;
	size_t len;
	struct nodelist *bqlist;
	int quotef;
	int oldstyle;
	/* syntax stack */
	struct synstack synbase = { .syntax = syntax };
	struct synstack *synstack = &synbase;

	if (syntax == DQSYNTAX)
		synstack->dblquote = 1;
	quotef = 0;
	bqlist = NULL;

	STARTSTACKSTR(out);
	loop: {	/* for each line, until end of word */
#if ATTY
		if (c == '\034' && doprompt
		 && attyset() && ! equal(termval(), "emacs")) {
			attyline();
			if (synstack->syntax == BASESYNTAX)
				return readtoken();
			c = pgetc_top(synstack);
			goto loop;
		}
#endif
		CHECKEND();	/* set c to PEOF if at end of here document */
		for (;;) {	/* until end of line or end of word */
			CHECKSTRSPACE(4, out);	/* permit 4 calls to USTPUTC */
			switch(synstack->syntax[c]) {
			case CNL:	/* '\n' */
				if (synstack->syntax == BASESYNTAX &&
				    !synstack->varnest)
					goto endword;	/* exit outer loop */
				USTPUTC(c, out);
				nlprompt();
				c = pgetc_top(synstack);
				goto loop;		/* continue outer loop */
			case CWORD:
				USTPUTC(c, out);
				break;
			case CCTL:
				if ((!eofmark) | synstack->dblquote |
				    synstack->varnest)
					USTPUTC(CTLESC, out);
				USTPUTC(c, out);
				break;
			/* backslash */
			case CBACK:
				c = pgetc2();
				if (c == PEOF) {
					USTPUTC(CTLESC, out);
					USTPUTC('\\', out);
					pungetc();
				} else {
					if (
						synstack->dblquote &&
						c != '\\' && c != '`' &&
						c != '$' && (
							c != '"' ||
							(eofmark != NULL &&
							 !synstack->varnest)
						) && (
							c != '}' ||
							!synstack->varnest
						)
					) {
						USTPUTC(CTLESC, out);
						USTPUTC('\\', out);
					}
					USTPUTC(CTLESC, out);
					USTPUTC(c, out);
					quotef++;
				}
				break;
			case CSQUOTE:
				synstack->syntax = SQSYNTAX;
quotemark:
				if (eofmark == NULL) {
					USTPUTC(CTLQUOTEMARK, out);
				}
				break;
			case CDQUOTE:
				synstack->syntax = DQSYNTAX;
				synstack->dblquote = 1;
toggledq:
				if (synstack->varnest)
					synstack->innerdq ^= 1;
				goto quotemark;
			case CENDQUOTE:
				if (eofmark && !synstack->varnest) {
					USTPUTC(c, out);
					break;
				}

				if (synstack->dqvarnest == 0) {
					synstack->syntax = BASESYNTAX;
					synstack->dblquote = 0;
				}

				quotef++;

				if (c == '"')
					goto toggledq;

				goto quotemark;
			case CVAR:	/* '$' */
				PARSESUB();		/* parse substitution */
				break;
			case CENDVAR:	/* '}' */
				if (!synstack->innerdq &&
				    synstack->varnest > 0) {
					if (!--synstack->varnest &&
					    synstack->varpushed)
						synstack_pop(&synstack);
					else if (synstack->dqvarnest > 0)
						synstack->dqvarnest--;
					USTPUTC(CTLENDVAR, out);
				} else {
					USTPUTC(c, out);
				}
				break;
			case CLP:	/* '(' in arithmetic */
				synstack->parenlevel++;
				USTPUTC(c, out);
				break;
			case CRP:	/* ')' in arithmetic */
				if (synstack->parenlevel > 0) {
					USTPUTC(c, out);
					--synstack->parenlevel;
				} else {
					if (pgetc_eatbnl() == ')') {
						USTPUTC(CTLENDARI, out);
						synstack_pop(&synstack);
					} else {
						/*
						 * unbalanced parens
						 *  (don't 2nd guess - no error)
						 */
						pungetc();
						USTPUTC(')', out);
					}
				}
				break;
			case CBQUOTE:	/* '`' */
				if (checkkwd & CHKEOFMARK) {
					USTPUTC('`', out);
					break;
				}

				PARSEBACKQOLD();
				break;
			case CEOF:
				goto endword;		/* exit outer loop */
			case CIGN:
				break;
			default:
				if (synstack->varnest == 0)
					goto endword;	/* exit outer loop */
				if (c != PEOA) {
					USTPUTC(c, out);
				}
			}
			c = pgetc_top(synstack);
		}
	}
endword:
	if (synstack->syntax == ARISYNTAX)
		synerror("Missing '))'");
	if (synstack->syntax != BASESYNTAX && eofmark == NULL)
		synerror("Unterminated quoted string");
	if (synstack->varnest != 0) {
		/* { */
		synerror("Missing '}'");
	}
	USTPUTC('\0', out);
	len = out - (char *)stackblock();
	out = stackblock();
	if (eofmark == NULL) {
		if ((c == '>' || c == '<')
		 && quotef == 0
		 && len <= 2
		 && (*out == '\0' || is_digit(*out))) {
			PARSEREDIR();
			return lasttoken = TREDIR;
		} else {
			pungetc();
		}
	}
	quoteflag = quotef;
	backquotelist = bqlist;
	grabstackblock(len);
	wordtext = out;
	return lasttoken = TWORD;
/* end of readtoken routine */



/*
 * Check to see whether we are at the end of the here document.  When this
 * is called, c is set to the first character of the next input line.  If
 * we are at the end of the here document, this routine sets the c to PEOF.
 */

checkend: {
	if (realeofmark(eofmark)) {
		int markloc;
		char *p;

		if (c == PEOA) {
			c = pgetc2();
		}
		if (striptabs) {
			while (c == '\t') {
				c = pgetc2();
			}
		}

		markloc = out - (char *)stackblock();
		for (p = eofmark; STPUTC(c, out), *p; p++) {
			if (c != *p)
				goto more_heredoc;

			c = pgetc2();
		}

		if (c == '\n' || c == PEOF) {
			c = PEOF;
			nlnoprompt();
		} else {
			int len;

more_heredoc:
			p = (char *)stackblock() + markloc + 1;
			len = out - p;

			if (len) {
				len -= c < 0;
				c = p[-1];

				if (len) {
					char *str;

					str = alloca(len + 1);
					*(char *)mempcpy(str, p, len) = 0;

					pushstring(str, NULL);
				}
			}
		}

		STADJUST((char *)stackblock() + markloc - out, out);
	}
	goto checkend_return;
}


/*
 * Parse a redirection operator.  The variable "out" points to a string
 * specifying the fd to be redirected.  The variable "c" contains the
 * first character of the redirection operator.
 */

parseredir: {
	char fd = *out;
	union node *np;

	np = (union node *)stalloc(sizeof (struct nfile));
	if (c == '>') {
		np->nfile.fd = 1;
		c = pgetc_eatbnl();
		if (c == '>')
			np->type = NAPPEND;
		else if (c == '|')
			np->type = NCLOBBER;
		else if (c == '&')
			np->type = NTOFD;
		else {
			np->type = NTO;
			pungetc();
		}
	} else {	/* c == '<' */
		np->nfile.fd = 0;
		switch (c = pgetc_eatbnl()) {
		case '<':
			if (sizeof (struct nfile) != sizeof (struct nhere)) {
				np = (union node *)stalloc(sizeof (struct nhere));
				np->nfile.fd = 0;
			}
			np->type = NHERE;
			heredoc = (struct heredoc *)stalloc(sizeof (struct heredoc));
			heredoc->here = np;
			if ((c = pgetc_eatbnl()) == '-') {
				heredoc->striptabs = 1;
			} else {
				heredoc->striptabs = 0;
				pungetc();
			}
			break;

		case '&':
			np->type = NFROMFD;
			break;

		case '>':
			np->type = NFROMTO;
			break;

		default:
			np->type = NFROM;
			pungetc();
			break;
		}
	}
	if (fd != '\0')
		np->nfile.fd = digit_val(fd);
	redirnode = np;
	goto parseredir_return;
}


/*
 * Parse a substitution.  At this point, we have read the dollar sign
 * and nothing else.
 */

parsesub: {
	int subtype;
	int typeloc;
	char *p;
	static const char types[] = "}-+?=";

	c = pgetc_eatbnl();
	if (
		(checkkwd & CHKEOFMARK) ||
		c <= PEOA  ||
		(c != '(' && c != '{' && !is_name(c) && !is_special(c))
	) {
		USTPUTC('$', out);
		pungetc();
	} else if (c == '(') {	/* $(command) or $((arith)) */
		if (pgetc_eatbnl() == '(') {
			PARSEARITH();
		} else {
			pungetc();
			PARSEBACKQNEW();
		}
	} else {
		const char *newsyn = synstack->syntax;

		USTPUTC(CTLVAR, out);
		typeloc = out - (char *)stackblock();
		STADJUST(1, out);
		subtype = VSNORMAL;
		if (likely(c == '{')) {
			c = pgetc_eatbnl();
			subtype = 0;
		}
varname:
		if (is_name(c)) {
			do {
				STPUTC(c, out);
				c = pgetc_eatbnl();
			} while (is_in_name(c));
		} else if (is_digit(c)) {
			do {
				STPUTC(c, out);
				c = pgetc_eatbnl();
			} while (is_digit(c));
		} else if (c != '}') {
			int cc = c;

			c = pgetc_eatbnl();

			if (!subtype && cc == '#') {
				subtype = VSLENGTH;

				if (c == '_' || isalnum(c))
					goto varname;

				cc = c;
				c = pgetc_eatbnl();
				if (cc == '}' || c != '}') {
					pungetc();
					subtype = 0;
					c = cc;
					cc = '#';
				}
			}

			if (!is_special(cc)) {
				if (subtype == VSLENGTH)
					subtype = 0;
				goto badsub;
			}

			USTPUTC(cc, out);
		} else
			goto badsub;

		if (subtype == 0) {
			int cc = c;

			switch (c) {
			case ':':
				subtype = VSNUL;
				c = pgetc_eatbnl();
				/*FALLTHROUGH*/
			default:
				p = strchr(types, c);
				if (p == NULL)
					break;
				subtype |= p - types + VSNORMAL;
				break;
			case '%':
			case '#':
				subtype = c == '#' ? VSTRIMLEFT :
						     VSTRIMRIGHT;
				c = pgetc_eatbnl();
				if (c == cc)
					subtype++;
				else
					pungetc();

				newsyn = BASESYNTAX;
				break;
			}
		} else {
badsub:
			pungetc();
		}

		if (newsyn == ARISYNTAX)
			newsyn = DQSYNTAX;

		if ((newsyn != synstack->syntax || synstack->innerdq) &&
		    subtype != VSNORMAL) {
			synstack_push(&synstack,
				      synstack->prev ?:
				      alloca(sizeof(*synstack)),
				      newsyn);

			synstack->varpushed++;
			synstack->dblquote = newsyn != BASESYNTAX;
		}

		*((char *)stackblock() + typeloc) = subtype;
		if (subtype != VSNORMAL) {
			synstack->varnest++;
			if (synstack->dblquote)
				synstack->dqvarnest++;
		}
		STPUTC('=', out);
	}
	goto parsesub_return;
}


/*
 * Called to parse command substitutions.  Newstyle is set if the command
 * is enclosed inside $(...); nlpp is a pointer to the head of the linked
 * list of commands (passed by reference), and savelen is the number of
 * characters on the top of the stack which must be preserved.
 */

parsebackq: {
	struct nodelist **nlpp;
	union node *n;
	char *str;
	size_t savelen;
	struct heredoc *saveheredoclist;
	int uninitialized_var(saveprompt);

	str = NULL;
	savelen = out - (char *)stackblock();
	if (savelen > 0) {
		str = alloca(savelen);
		memcpy(str, stackblock(), savelen);
	}
        if (oldstyle) {
                /* We must read until the closing backquote, giving special
                   treatment to some slashes, and then push the string and
                   reread it as input, interpreting it normally.  */
                char *pout;
                int pc;
                size_t psavelen;
                char *pstr;


                STARTSTACKSTR(pout);
		for (;;) {
			if (needprompt) {
				setprompt(2);
			}
			switch (pc = pgetc_eatbnl()) {
			case '`':
				goto done;

			case '\\':
                                pc = pgetc_eatbnl();
                                if (pc != '\\' && pc != '`' && pc != '$'
                                    && (!synstack->dblquote || pc != '"'))
                                        STPUTC('\\', pout);
				if (pc > PEOA) {
					break;
				}
				/* fall through */

			case PEOF:
			case PEOA:
				synerror("EOF in backquote substitution");

			case '\n':
				nlnoprompt();
				break;

			default:
				break;
			}
			STPUTC(pc, pout);
                }
done:
                STPUTC('\0', pout);
                psavelen = pout - (char *)stackblock();
                if (psavelen > 0) {
			pstr = grabstackstr(pout);
			setinputstring(pstr);
                }
        }
	nlpp = &bqlist;
	while (*nlpp)
		nlpp = &(*nlpp)->next;
	*nlpp = (struct nodelist *)stalloc(sizeof (struct nodelist));
	(*nlpp)->next = NULL;

	saveheredoclist = heredoclist;
	heredoclist = NULL;

	if (oldstyle) {
		saveprompt = doprompt;
		doprompt = 0;
	}

	n = list(2);

	if (oldstyle)
		doprompt = saveprompt;
	else {
		if (readtoken() != TRP)
			synexpect(TRP);
		setinputstring(nullstr);
		parseheredoc();
	}

	heredoclist = saveheredoclist;

	(*nlpp)->n = n;
	/* Start reading from old file again. */
	popfile();
	/* Ignore any pushed back tokens left from the backquote parsing. */
	if (oldstyle)
		tokpushback = 0;
	out = growstackto(savelen + 1);
	if (str) {
		memcpy(out, str, savelen);
		STADJUST(savelen, out);
	}
	USTPUTC(CTLBACKQ, out);
	if (oldstyle)
		goto parsebackq_oldreturn;
	else
		goto parsebackq_newreturn;
}

/*
 * Parse an arithmetic expansion (indicate start of one and set state)
 */
parsearith: {

	synstack_push(&synstack,
		      synstack->prev ?: alloca(sizeof(*synstack)),
		      ARISYNTAX);
	synstack->dblquote = 1;
	USTPUTC(CTLARI, out);
	goto parsearith_return;
}

} /* end of readtoken */



#ifdef mkinit
INCLUDE "parser.h"
#endif


/*
 * Return of a legal variable name (a letter or underscore followed by zero or
 * more letters, underscores, and digits).
 */

char *
endofname(const char *name)
	{
	char *p;

	p = (char *) name;
	if (! is_name(*p))
		return p;
	while (*++p) {
		if (! is_in_name(*p))
			break;
	}
	return p;
}


/*
 * Called when an unexpected token is read during the parse.  The argument
 * is the token that is expected, or -1 if more than one type of token can
 * occur at this point.
 */

STATIC void
synexpect(int token)
{
	char msg[64];

	if (token >= 0) {
		fmtstr(msg, 64, "%s unexpected (expecting %s)",
			tokname[lasttoken], tokname[token]);
	} else {
		fmtstr(msg, 64, "%s unexpected", tokname[lasttoken]);
	}
	synerror(msg);
	/* NOTREACHED */
}


STATIC void
synerror(const char *msg)
{
	errlinno = plinno;
	sh_error("Syntax error: %s", msg);
	/* NOTREACHED */
}

STATIC void
setprompt(int which)
{
	struct stackmark smark;
	int show;

	needprompt = 0;
	whichprompt = which;

#ifdef SMALL
	show = 1;
#else
	show = !el;
#endif
	if (show) {
		pushstackmark(&smark, stackblocksize());
		out2str(getprompt(NULL));
		popstackmark(&smark);
	}
}

const char *
expandstr(const char *ps)
{
	union node n;
	int saveprompt;

	/* XXX Fix (char *) cast. */
	setinputstring((char *)ps);

	saveprompt = doprompt;
	doprompt = 0;

	readtoken1(pgetc_eatbnl(), DQSYNTAX, FAKEEOFMARK, 0);

	doprompt = saveprompt;

	popfile();

	n.narg.type = NARG;
	n.narg.next = NULL;
	n.narg.text = wordtext;
	n.narg.backquote = backquotelist;

	expandarg(&n, NULL, EXP_QUOTED);
	return stackblock();
}

/*
 * called by editline -- any expansions to the prompt
 *    should be added here.
 */
const char *
getprompt(void *unused)
{
	const char *prompt;

	switch (whichprompt) {
	default:
#ifdef DEBUG
		return "<internal prompt error>";
#endif
	case 0:
		return nullstr;
	case 1:
		prompt = ps1val();
		break;
	case 2:
		prompt = ps2val();
		break;
	}

	return expandstr(prompt);
}

const char *const *
findkwd(const char *s)
{
	return findstring(
		s, parsekwd, sizeof(parsekwd) / sizeof(const char *)
	);
}
