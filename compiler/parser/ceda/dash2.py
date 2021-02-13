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
