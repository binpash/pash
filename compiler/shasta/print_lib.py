STRING_OF_VAR_TYPE_DICT = {
    "Normal"   : "",
    "Minus"    : "-",
    "Plus"     : "+",
    "Question" : "?",
    "Assign"   : "=",
    "TrimR"    : "%",
    "TrimRMax" : "%%",
    "TrimL"    : "#",
    "TrimLMax" : "##",
    "Length"   : "#"
}

UNQUOTED = 0 # everything escaped
QUOTED = 1   # only escape special characters
HEREDOC = 2  # like QUOTED, but _don't_ escape double quotes
QUOTE_MODES = [UNQUOTED, QUOTED, HEREDOC]

def intercalate(p, ss):
    str = p.join(ss)
    return(str)

def braces(s):
    return "{ " + s + " ; }"

def parens(s):
    return f"( {s} )"

def string_of_var_type(var_type):
    if var_type in STRING_OF_VAR_TYPE_DICT:
        return STRING_OF_VAR_TYPE_DICT[var_type]
    exit(1)

## TODO: Fix this
def separated(f, l):
    return " ".join(map(f, l))

def show_unless(expected, actual):
    if expected == actual:
        return ""
    else:
        return str(actual)

def background(s):
    return "{ " + s + " & }"

def string_of_redirs(rs):
    str = ""

    for redir in rs:
        str = str + " " + redir.pretty()

    return str

def fresh_marker (heredoc):
    respectsFound = set();

    for line in heredoc.split ('\n'):
        respects = 0;

        if ((len (line) > 2) and (line [0] == 'E') and (line [1] == 'O')):
            for i in range (2, len (line)):
                if (line [i] == 'F'):
                    respects = i - 2;

            respectsFound.add(respects);

    i = 0;
    while (True):
        if (not (i in respectsFound)):
            return "EOF" + ("F" * i);

        i = i + 1;


# This version may give an unnecessarily long EOFFFF... (and therefore won't
# match the OCaml output but it is still correct w.r.t. giving a fresh
# marker, and uses less memory than fresh_marker above).
def fresh_marker0(heredoc):
    maxRespects = 0

    for line in heredoc.split('\n'):
        respects = 0

        if (len (line) > 2) and (line [0] == 'E') and (line [1] == 'O'):
            for i in range(2, len (line)):
                if line[i] == 'F':
                    respects = i - 1

            maxRespects = max(maxRespects, respects)

    return "EOF" + ("F" * maxRespects)

def escaped(param):
    char = chr(param)

    if char == "'":
        return "\\'"
    elif char == "\\":
        return "\\\\"
    elif char == "\n":
        return "\\n"
    elif char == "\t":
        return "\\t"
    elif char == "\r":
        return "\\r"
    elif char == "\b":
        return "\\b"
    elif (param >= ord (' ')) and (param <= ord ('~')):
        return char
    else:
#        str1 =   "\\" \
#               + chr (48 + int (param / 100)) \
#               + chr (48 + ((int (param / 10)) % 10)) \
#               + chr (48 + (param % 10));
        return "\\" + str (param)

