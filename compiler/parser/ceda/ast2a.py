# import os;
import sys;
from dash2 import *;


# parser.h
CTLESC       = 129;
CTLVAR       = 130;
CTLENDVAR    = 131;
CTLBACKQ     = 132;
CTLARI       = 134;
CTLENDARI    = 135;
CTLQUOTEMARK = 136;

# Internal use only
STACK_CTLVAR = 100;
STACK_CTLARI = 101;
STACK_CTLQUO = 102;


VAR_TYPES \
    = [
       "Normal",   # 0x0
       "UNUSED",
       "Minus",    # 0x2
       "Plus",     # 0x3
       "Question", # 0x4
       "Assign",   # 0x5
       "TrimR",    # 0x6
       "TrimRMax", # 0x7
       "TrimL",    # 0x8
       "TrimLMax", # 0x9
       "Length"    # 0xa
      ];


SKIP_COMMAND = ["Command", [-1, [], [], []]];


def var_type (i):
    return VAR_TYPES [i];


def of_node (n_ptr):
    if (not n_ptr):
        return SKIP_COMMAND;
    elif (n_ptr == None):
        print ("Didn't expect this type of null");
        os.abort ();
    else:
        n = n_ptr.contents;

        if (n.type == NCMD):
            return (["Command",
                     [n.ncmd.linno,
                      to_assigns (n.ncmd.assign),
                      to_args (n.ncmd.args),
                      redirs (n.ncmd.redirect)]]);
        elif (n.type == NPIPE):
            return (["Pipe",
                     [n.npipe.backgnd != 0,
                      list (map (of_node, nodelist (n.npipe.cmdlist)))]]);
        elif (n.type == NREDIR):
            return ["Redir", of_nredir (n)];
        elif (n.type == NBACKGND):
            return ["Background", of_nredir (n)];
        elif (n.type == NSUBSHELL):
            return ["Subshell", of_nredir (n)];
        elif (n.type == NAND):
            return ["And", of_binary (n)];
        elif (n.type == NOR):
            return ["Or", of_binary (n)];
        elif (n.type == NSEMI):
            return ["Semi", of_binary (n)];
        elif (n.type == NIF):
            return (["If",
                     [of_node (n.nif.test),
                      of_node (n.nif.ifpart),
                      of_node (n.nif.elsepart)]]);
        elif (n.type == NWHILE):
            return ["While", of_binary (n)];
        elif (n.type == NUNTIL):
            (t, b) = of_binary (n);
            return ["While", [["Not", t], b]];
        elif (n.type == NFOR):
            return ["For",
                    [n.nfor.linno,
                     to_arg (n.nfor.args.contents.narg),
                     of_node (n.nfor.body),
                     n.nfor.var.decode ("charmap")]];
        elif (n.type == NCASE):
            cases_hashes = []; # Poetic

            for case in caselist (n.ncase.cases):
                (pattern, body) = case;

                current_case \
                    = {'cpattern' : to_args (pattern),
                       'cbody'    : of_node (body)};

                cases_hashes.append (current_case);

            return (["Case",
                     [n.ncase.linno,
                      to_arg (n.ncase.expr.contents.narg),
                      cases_hashes]]);
        elif (n.type == NDEFUN):
            return ["Defun",
                    [n.ndefun.linno,
                     n.ndefun.text.decode ("charmap"),
                     of_node (n.ndefun.body)]];
        elif (n.type == NNOT):
            return ["Not", of_node (n.nnot.com)];
        else:
            print ("Unexpected type");
            sys.stdout.flush ();
            os.abort ();


def of_nredir (n):
    return ([n.nredir.linno, of_node (n.nredir.n), redirs (n.nredir.redirect)]);


def mk_file (ty, n):
    arg = to_arg (n.nfile.fname.contents.narg);

    return ["File", [ty, n.nfile.fd, arg]];


def mk_dup (ty, n):
    ndup = n.ndup;
    vname = ndup.vname;

    tgt = [];

    if (not vname):
        dupfd = ndup.dupfd;
        if (dupfd == -1):
            tgt.append (["C", ord ("-")]);
        else:
            dupfd_str = str (dupfd);

            for i in range (len (dupfd_str)):
                tgt.append (["C", ord (dupfd_str [i])]);
    else:
        tgt = to_arg (vname.narg);

    return (["Dup", [ty, ndup.fd, tgt]]);


def mk_here (ty, n):
    return ["Heredoc", [ty, n.nhere.fd, to_arg (n.nhere.doc.contents.narg)]];


def redirs (n_ptr):
    rlist = [];

    while (n_ptr):
        h = [];

        n = n_ptr.contents;

        if (n.type == NTO):
            h = mk_file ("To", n);
        elif (n.type == NCLOBBER):
            h = mk_file ("Clobber", n);
        elif (n.type == NFROM):
            h = mk_file ("From", n);
        elif (n.type == NFROMTO):
            h = mk_file ("FromTo", n);
        elif (n.type == NAPPEND):
            h = mk_file ("Append", n);
        elif (n.type == NTOFD):
            h = mk_dup ("ToFD", n);
        elif (n.type == NFROMFD):
            h = mk_dup ("FromFD", n);
        elif (n.type == NHERE):
            h = mk_here ("Here", n);
        elif (n.type == NXHERE):
            h = mk_here ("XHere", n);
        else:
            print ("unexpected node_type in redirlist");
            os.abort ();

        rlist.append (h);

        n_ptr = n.nfile.next;

    return rlist;


def of_binary (n):
    return [of_node (n.nbinary.ch1), of_node (n.nbinary.ch2)];


def to_arg (narg):
    arg = to_arg_rev (narg);
    arg.reverse ();

    return (arg);


def to_arg_rev (narg):
    (a, s, bqlist, stack) = parse_arg (explode_rev (narg.text), narg.backquote, []);

    assert (len (s) == 0);
    # assert (nullptr bqlist)
    if (bqlist):
        print ("bqlist is not null");
        print (bqlist);
        os.abort ();
    assert (len (stack) == 0);

    return (a);


def parse_arg (s, bqlist, stack):
    # | [],[] -> [],[],bqlist,[]
    if ((len (s) == 0) and (len (stack) == 0)):
        return ([], [], bqlist, []);
    # | [],`CTLVar::_ -> failwith "End of string before CTLENDVAR"
    elif ((len (s) == 0) and (len (stack) > 0) and (stack [-1] == STACK_CTLVAR)):
        print ("End of string before CTLENDVAR");
        os.abort ();
    # | [],`CTLAri::_ -> failwith "End of string before CTLENDARI"
    elif ((len (s) == 0) and (len (stack) > 0) and (stack [-1] == STACK_CTLARI)):
        print (s);
        print (stack);

        print ("End of string before CTLENDARI");
        os.abort ();
    # | [],`CTLQuo::_ -> failwith "End of string before CTLQUOTEMARK"
    elif ((len (s) == 0) and (len (stack) > 0) and (stack [-1] == STACK_CTLQUO)):
        print (s);
        print (stack);

        print ("End of string before CTLENDQUOTEMARK");
        os.abort ();

    # (* CTLESC *)
    # | '\129'::c::s,_ -> arg_char (E c) s bqlist stack
    elif ((len (s) >= 2) and (s [-1] == CTLESC)):
        s.pop ();
        c = s.pop ();

        return (arg_char (["E", c], s, bqlist, stack));

    # (* CTLVAR *)
    # | '\130'::t::s,_ ->
    elif ((len (s) >= 2) and (s [-1] == CTLVAR)):
        s.pop ();
        t = s.pop ();

        # let var_name,s = split_at (fun c -> c = '=') s in
        var_name = []
        while ((len (s) > 0) and (s [-1] != ord ('='))):
            c = s.pop ();
            var_name.append (c);

        v = [];

        if (((t & 0xf) == 0x1) and (len (s) >= 1) and (s [-1] == ord ('='))):
            s.pop ();

            v = ["V", ["Normal", False, implode (var_name), []]];
        elif (((t & 0xf) == 0xa) and (len (s) >= 2) and (s [-1] == ord ('=')) and (s [-2] == 131)):
            s.pop ();
            s.pop ();

            v = ["V", ["Length", False, implode (var_name), []]];
        elif ((((t & 0xf) == 0x1) or ((t & 0xf) == 0xa)) and (len (s) >= 1)):
            print ("Missing CTLENDVAR for VSNORMAL/VSLENGTH");
            os.abort ();
        elif ((len (s) >= 1) and (s [-1] == ord ('='))):
            s.pop ();

            vstype = t & 0xf;

            stack.append (STACK_CTLVAR);

            (a, s, bqlist, stack) = parse_arg (s, bqlist, stack);
            a.reverse ();

            v = ["V", [var_type (vstype), (t & 0x10 == 0x10), implode (var_name), a]];
        elif (len (s) >= 1):
            print (s);
            print (stack);

            print ("Expected '=' terminating variable name");
            os.abort ();
        elif (len (s) == 0):
            print ("Expected '=' terminating variable name, found EOF");
            os.abort ();
        else:
            print ("This shouldn't be reachable");
            os.abort ();

        return arg_char (v, s, bqlist, stack);
    # | '\130'::_, _ -> raise (ParseException "bad substitution (missing variable name in ${}?")
    elif (False and (len (s) >= 1) and (s [-1] == CTLVAR)):
        print (s);
        print (stack);

        print ("bad substitution (missing variable name in ${}?");
        os.abort ();

    # (* CTLENDVAR *)
    # | '\131'::s,`CTLVar::stack' -> [],s,bqlist,stack'
    elif ((len (s) >= 1) and (s [-1] == CTLENDVAR) and (len (stack) >= 1) and (stack [-1] == STACK_CTLVAR)):
        s.pop ();
        stack.pop ();

        return ([], s, bqlist, stack);
    # | '\131'::_,`CTLAri::_ -> failwith "Saw CTLENDVAR before CTLENDARI"
    elif ((len (s) >= 1) and (s [-1] == CTLENDVAR) and (len (stack) >= 1) and (stack [-1] == STACK_CTLARI)):
        print ("Saw CTLENDVAR before CTLENDARI");
        os.abort ();
    # | '\131'::_,`CTLQuo::_ -> failwith "Saw CTLENDVAR before CTLQUOTEMARK"
    elif ((len (s) >= 1) and (s [-1] == CTLENDVAR) and (len (stack) >= 1) and (stack [-1] == STACK_CTLQUO)):
        print ("Saw CTLENDVAR before CTLQUOTEMARK");
        os.abort ();
    # | '\131'::_,[] -> failwith "Saw CTLENDVAR outside of CTLVAR"
    elif ((len (s) >= 1) and (s [-1] == CTLENDVAR)):
        print ("Saw CTLENDVAR outside of CTLVAR");
        os.abort ();

    # (* CTLBACKQ *)
    # | '\132'::s,_ ->
    elif ((len (s) >= 1) and (s [-1] == CTLBACKQ)):
        s.pop ();

        if (not bqlist):
            print (bqlist);
            print ("Saw CTLBACKQ but bqlist was null");
            os.abort ();
        else:
            return arg_char (["B", of_node (bqlist.contents.n)], s, bqlist.contents.next, stack);

    # (* CTLARI *)
    # | '\134'::s,_ ->
    elif ((len (s) >= 1) and (s [-1] == CTLARI)):
        s.pop ();

        stack.append (STACK_CTLARI);

        (a, s, bqlist, stack) = parse_arg (s, bqlist, stack);
        a.reverse ();

        # TODO: assert (stack = stack');

        return arg_char (["A", a], s, bqlist, stack);

    # (* CTLENDARI *)
    # | '\135'::s,`CTLAri::stack' -> [],s,bqlist,stack'
    elif ((len (s) >= 1) and (s [-1] == CTLENDARI) and (len (stack) >= 1) and (stack [-1] == STACK_CTLARI)):
        s.pop ();
        stack.pop ();

        return ([], s, bqlist, stack);
    # | '\135'::_,`CTLVar::_' -> failwith "Saw CTLENDARI before CTLENDVAR"
    elif ((len (s) >= 1) and (s [-1] == CTLENDARI) and (len (stack) >= 1) and (stack [-1] == STACK_CTLVAR)):
        print ("Saw CTLENDARI before CTLENDVAR");
        os.abort ();
    # | '\135'::_,`CTLQuo::_' -> failwith "Saw CTLENDARI before CTLQUOTEMARK"
    elif ((len (s) >= 1) and (s [-1] == CTLENDARI) and (len (stack) >= 1) and (stack [-1] == STACK_CTLQUO)):
        print ("Saw CTLENDARI before CTLQUOTEMARK");
        os.abort ();
    # | '\135'::_,[] -> failwith "Saw CTLENDARI outside of CTLARI"
    elif ((len (s) >= 1) and (s [-1] == CTLENDARI) and (len (stack) == 0)):
        print ("Saw CTLENDARI outside of CTLARI");
        os.abort ();

    # (* CTLQUOTEMARK *)
    # | '\136'::s,`CTLQuo::stack' -> [],s,bqlist,stack'
    elif ((len (s) >= 1) and (s [-1] == CTLQUOTEMARK) and (len (stack) >= 1) and (stack [-1] == STACK_CTLQUO)):
        s.pop ();
        stack.pop ();

        return ([], s, bqlist, stack);
    # | '\136'::s,_ ->
    elif ((len (s) >= 1) and (s [-1] == CTLQUOTEMARK)):
        s.pop ();
        stack.append (STACK_CTLQUO);

        (a, s, bqlist, stack) = parse_arg (s, bqlist, stack);
        a.reverse ();

        return arg_char (["Q", a], s, bqlist, stack);

    # (* tildes *)
    # | '~'::s,stack ->
    elif ((len (s) >= 1) and (s [-1] == ord ('~'))):
        s.pop ();

        if ((STACK_CTLQUO in stack) or (STACK_CTLARI in stack)):
            return arg_char (["C", "~"], s, bqlist, stack);
        else:
#            print ("Calling parse_tilde");
#            print (s);

            (uname, sL) = parse_tilde ([], s);

            return arg_char (["T", uname], sL, bqlist, stack);
    # (* ordinary character *)
    # | c::s,_ -> arg_char (C c) s bqlist stack
    elif (len (s) >= 1):
        c = s.pop ();

        return (arg_char (["C", c], s, bqlist, stack));
    else:
        print ("parse_arg: unreachable case");

        print (s);
        print (stack);

        os.abort ();


def implodeOrNull (acc):
    if (acc == []):
        return "None";
    else:
        acc.reverse ();
        return ["Some", implode (acc)];


def parse_tilde (acc, s):
    if (s == []):
#        print ("Acc 0: ");
        print (acc);
        return (implodeOrNull (acc), []);
    elif ((len (s) >= 1) and (s [-1] == 129)):
        return (["None"], s);
    elif ((len (s) >= 1) and (s [-1] == 136)):
        return (["None"], s);
    elif ((len (s) >= 1) and (s [-1] == 131)):
#        print ("Acc 1: ");
#        print (acc);

        return (implodeOrNull (acc), s);
    elif ((len (s) >= 1) and (s [-1] == ord (':'))):
#        print ("Acc 2: ");
#        print (acc);

        return (implodeOrNull (acc), s);
    elif ((len (s) >= 1) and (s [-1] == ord ('/'))):
#        print ("Acc 3: ");
#        print (acc);

        return (implodeOrNull (acc), s);
    else:
        c = s.pop ();
        acc.append (c);

        return parse_tilde (acc, s);


# arg_char c s bqlist stack =
#  let a,s,bqlist,stack = parse_arg s bqlist stack in
#  (c::a,s,bqlist,stack)
def arg_char (c, s, bqlist, stack):
    (a, s, bqlist, stack) = parse_arg (s, bqlist, stack);

    a.append (c);

    return (a, s, bqlist, stack);


def to_assign (v_rev, a_rev):
    if (len (a_rev) == 0):
        print ("Never found an '=' sign in assignment");
        os.abort ();
    elif ((a_rev [-1][0] == 'C') and (a_rev [-1][1] == ord ('='))):
        a_rev.pop ();

        a_rev.reverse ();
        return (implode (v_rev), a_rev);
    elif ((a_rev [-1][0] == 'C')):
        c = a_rev [-1][1];
        a_rev.pop ();
        v_rev.append (c);

        return (to_assign (v_rev, a_rev));
    else:
        print ("Unexpected special character in assignment");
        sys.stdout.flush ();
        os.abort ();


# to_assigns n = List.map (to_assign []) (to_args n)
def to_assigns (n):
    assigns = []

    for i in (to_args (n)):
        i.reverse ();
        assigns.append (to_assign ([], i));

    assigns.reverse ();
    return (assigns);


def to_args (n):
    snek = [];

    # ctypes has different semantics for POINTER vs. c_void_p
    # See https://groups.google.com/g/nzpug/c/5CJxaWjuQro
    while (n):
        if (n.contents.type != NARG):
            print ("Unexpected type: " + n.contents.type);
            sys.stdout.flush ();
            os.abort ();

        arg = to_arg (n.contents.narg);
        snek.append (arg);
        n = n.contents.narg.next;

    return snek;
