from ctypes import *


# nodes.h
NCMD      = 0;
NPIPE     = 1;
NREDIR    = 2;
NBACKGND  = 3;
NSUBSHELL = 4;
NAND      = 5;
NOR       = 6;
NSEMI     = 7;
NIF       = 8;
NWHILE    = 9;
NUNTIL    = 10;
NFOR      = 11;
NCASE     = 12;
NCLIST    = 13;
NDEFUN    = 14;
NARG      = 15;
NTO       = 16;
NCLOBBER  = 17;
NFROM     = 18;
NFROMTO   = 19;
NAPPEND   = 20;
NTOFD     = 21;
NFROMFD   = 22;
NHERE     = 23;
NXHERE    = 24;
NNOT      = 25;


# struct stackmark {
#     struct stack_block *stackp;
#     char *stacknxt;
#     size_t stacknleft;
# };
#
# We only care about getting the struct size correct, not the contents.
class stackmark (Structure):
    _fields_ = [("stackp", c_void_p),
                ("nxt",    c_void_p),
                ("size",   c_size_t)];

def init_stack (libdash):
    stack = create_string_buffer (sizeof (stackmark));

    libdash.setstackmark.argtypes = [c_void_p]; # Pretend we don't know the contents
    libdash.setstackmark.restypes = None;
    libdash.setstackmark (stack);

    return (stack);

def pop_stack (libdash, smark):
    # Inefficient, we should only initialize once

    libdash.popstackmark.argtypes = [c_void_p]; # Again, hide the contents
    libdash.popstackmark.restype = None;

    return (libdash.popstackmark (smark));


def dash_init (libdash):
    libdash.init.argtypes = [];
    libdash.init.restype = None;

    libdash.init ();


def initialize_dash_errno (libdash):
    libdash.initialize_dash_errno.argtypes = [];
    libdash.initialize_dash_errno.restype = None;

    libdash.initialize_dash_errno ();


def initialize (libdash):
    initialize_dash_errno (libdash);
    dash_init (libdash);


def setinputtostdin (libdash):
    libdash.setinputfd.argtypes = [c_int, c_int];
    libdash.setinputfd.restype = None;

    libdash.setinputfd (0, 0);


# TODO: allow push parameter
def setinputfile (libdash, filename):
    libdash.setinputfile.argtypes = [c_char_p, c_int];
    libdash.setinputfile.restypes = c_int;
    libdash.setinputfile (filename.encode ('utf-8'), 0);


def parsecmd_safe (libdash, interactive):
    libdash.parsecmd_safe.argtypes = [c_int];
    libdash.parsecmd_safe.restype = c_void_p;

    return (libdash.parsecmd_safe (int (interactive)));


# Forward declarations to break recursive dependencies
class union_node (Union):
    pass;

class nodelist (Structure):
    pass;


class ncmd (Structure):
    _fields_ = [("type",     c_int),
                ("linno",    c_int),
                ("assign",   POINTER (union_node)),
                ("args",     POINTER (union_node)),
                ("redirect", POINTER (union_node))];

class npipe (Structure):
    _fields_ = [("type",   c_int),
               ("backgnd", c_int),
               ("cmdlist", POINTER (nodelist))];

class nredir (Structure):
    _fields_ = [("type",     c_int),
                ("linno",    c_int),
                ("n",        POINTER (union_node)),
                ("redirect", POINTER (union_node))];

class nbinary (Structure):
    _fields_ = [("type", c_int),
                ("ch1",  POINTER (union_node)),
                ("ch2",  POINTER (union_node))];

class nif (Structure):
    _fields_ = [("type",     c_int),
                ("test",     POINTER (union_node)),
                ("ifpart",   POINTER (union_node)),
                ("elsepart", POINTER (union_node))];

class nfor (Structure):
    _fields_ = [("type",  c_int),
                ("linno", c_int),
                ("args",  POINTER (union_node)),
                ("body",  POINTER (union_node)),
                ("var",   c_char_p)];

class ncase (Structure):
    _fields_ = [("type",  c_int),
                ("linno", c_int),
                ("expr",  POINTER (union_node)),
                ("cases", POINTER (union_node))];

class nclist (Structure):
    _fields_ = [("type",    c_int),
                ("next",    POINTER (union_node)),
                ("pattern", POINTER (union_node)),
                ("body",    POINTER (union_node))];

class ndefun (Structure):
    _fields_ = [("type",  c_int),
                ("linno", c_int),
                ("text",  c_char_p),
                ("body",  POINTER (union_node))];

class narg (Structure):
    _fields_ = [("type",      c_int),
                ("next",      POINTER (union_node)),
                ("text",      c_char_p),
                ("backquote", POINTER (nodelist))];

class nfile (Structure):
    _fields_ = [("type",     c_int),
                ("next",     POINTER (union_node)),
                ("fd",       c_int),
                ("fname",    POINTER (union_node)),
                ("expfname", c_char_p)]

class ndup (Structure):
    _fields_ = [("type",  c_int),
                ("next",  POINTER (union_node)),
                ("fd",    c_int),
                ("dupfd", c_int),
                ("vname", POINTER (union_node))];

class nhere (Structure):
    _fields_ = [("type", c_int),
                ("next", POINTER (union_node)),
                ("fd",   c_int),
                ("doc",  POINTER (union_node))];

class nnot (Structure):
    _fields_ = [("type", c_int),
                ("com",  POINTER (union_node))];


nodelist._fields_ = [("next", POINTER (nodelist)),
                     ("n",    POINTER (union_node))];

union_node._fields_ = [("type",    c_int),
                       ("ncmd",    ncmd),
                       ("npipe",   npipe),
                       ("nredir",  nredir),
                       ("nbinary", nbinary),
                       ("nif",     nif),
                       ("nfor",    nfor),
                       ("ncase",   ncase),
                       ("nclist",  nclist),
                       ("ndefun",  ndefun),
                       ("narg",    narg),
                       ("nfile",   nfile),
                       ("ndup",    ndup),
                       ("nhere",   nhere),
                       ("nnot",    nnot)];


class strpush (Structure):
    pass;

# struct strpush {
#     struct strpush *prev;   /* preceding string on stack */
#     char *prevstring;
#     int prevnleft;
#     struct alias *ap;   /* if push was associated with an alias */
#     char *string;       /* remember the string since it may change */
#
#     /* Remember last two characters for pungetc. */
#     int lastc[2];
#
#     /* Number of outstanding calls to pungetc. */
#     int unget;
# };
strpush._fields_ = [("prev", POINTER (strpush)),
                    ("prevstring", c_char_p),
                    ("prevnleft", c_int),
                    ("ap", c_void_p),
                    ("string", c_char_p),
                    ("lastc", 2 * c_int),
                    ("unget", c_int)];

class parsefile (Structure):
    pass;

# struct parsefile {
#     struct parsefile *prev; /* preceding file on stack */
#     int linno;      /* current line */
#     int fd;         /* file descriptor (or -1 if string) */
#     int nleft;      /* number of chars left in this line */
#     int lleft;      /* number of chars left in this buffer */
#     char *nextc;        /* next char in buffer */
#     char *buf;      /* input buffer */
#     struct strpush *strpush; /* for pushing strings at this level */
#     struct strpush basestrpush; /* so pushing one is fast */
#
#     /* Remember last two characters for pungetc. */
#     int lastc[2];
#
#     /* Number of outstanding calls to pungetc. */
#     int unget;
# };
parsefile._fields_ = [("prev",        POINTER (parsefile)),
                      ("linno",       c_int),
                      ("fd",          c_int),
                      ("nleft",       c_int),
                      ("lleft",       c_int),
                      ("nextc",       POINTER (c_char)), # NOT c_char_p!
                      ("buf",         c_char_p),
                      ("strpush",     POINTER (strpush)),
                      ("basestrpush", strpush),
                      ("lastc",       2 * c_int),
                      ("unget",       c_int)];


# dash.ast
# let rec nodelist (n : nodelist structure ptr) : (node union ptr) list =
#   if nullptr n
#   then []
#   else (n @-> nodelist_n)::nodelist (n @-> nodelist_next)
def nodelist (nl):
    snek = [];

    # ctypes has different semantics for POINTER vs. c_void_p
    # See https://groups.google.com/g/nzpug/c/5CJxaWjuQro
    while (nl):
        snek.append (nl.contents.n);
        nl = nl.contents.next;

    return snek;


def caselist (n):
    cases = [];

    while (n):
        nclist = n.contents.nclist;

        assert (nclist.type == 13);

        cases.append ((nclist.pattern, nclist.body));

        n = nclist.next;

    return (cases);


def explode_rev (bytes):
    charlist = explode (bytes);
    charlist.reverse ();

    return (charlist);


def explode (bytes):
    s = bytes.decode ("charmap");

    charlist = [];

    for i in range (len (s)):
        charlist.append (ord (s [i]));

    return (charlist);


def implode_rev (l):
    s = implode (reversed (l));

    return (s);


def implode (l):
    s = "";

    for c in l:
        s = s + chr (c);

    return (s);
