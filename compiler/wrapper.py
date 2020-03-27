orig_commands = []
cid = 0
WRAPPER_CMD = "wrap"

def get_results():
    global orig_commands, cid, WRAPPER_CMD
    return (cid, orig_commands, WRAPPER_CMD)
    
def get_id():
    global cid 
    sid = str(cid)
    cid = cid + 1
    return sid

def get_wrapper():
    # ord / chr() for ascii conversion
    global WRAPPER_CMD
    io = [ ['C', ord(c)] for c in WRAPPER_CMD]
    id = [ ['C', ord(c)] for c in get_id()]
    return (io, id)

def wrap_command_value(value):
    global orig_commands
    orig_commands.append(value)
    w = get_wrapper()
    value[2].insert(0, w[1])
    value[2].insert(0, w[0])

def try_wrap(ast_node):
    # FIXME: deprecated!
    if ast_node[0] == "Command":
        value = ast_node[1]
        # For commands, `value` is of type:
        # [ lineno, [assgn], [arg], [redir] ]
        wrap_command_value(value)

def rewrite_ast(asts):
    for i, node in enumerate(asts):
        # print("Rewrite AST {}".format(i))
        rewrite_node(node)

# Fixme: Pythonic conversion for in-place rewrite? e.g., `rewrite!` ?
def rewrite_node(node):
    rewrite = {
        "Command": rewrite_command,
        "Pipe": rewrite_pipe,
        "Redir": rewrite_redir,
        "Background": rewrite_background,
        "Subshell": rewrite_subshell,
        "And": rewrite_and,
        "Or": rewrite_or,
        "Not": rewrite_not,
        "Semi": rewrite_semi,
        "If": rewrite_if,
        "While": rewrite_while,
        "For": rewrite_for,
        "Case": rewrite_case,
        "Defun": rewrite_defun
    }
    rewrite[node[0]](node[1])

def rewrite_command(value):
    # Command of (linno * assign list * args * redirection list)
    wrap_command_value(value)

def rewrite_pipe(value):
    # Pipe of (bool * t list)
    for i, term in enumerate(value[1]):
        rewrite_node(term)
    
def rewrite_redir(value):
    # Redir of (linno * t * redirection list)
    rewrite_node(value[1])

def rewrite_background(value):
    # Background of (linno * t * redirection list)
    rewrite_node(value[1])

def rewrite_subshell(value):
    # Subshell of (linno * t * redirection list)
    rewrite_node(value[1])

def rewrite_and(value):
    # And of (t * t)
    rewrite_node(value[0])
    rewrite_node(value[1])

def rewrite_or(value):
    # Or of (t * t)
    rewrite_node(value[0])
    rewrite_node(value[1])

def rewrite_not(value):
    # Not of t
    rewrite_node(value[0])

def rewrite_semi(value):
    # Semi of (t * t)
    rewrite_node(value[0])
    rewrite_node(value[1])

def rewrite_if(value):
    # If of (t * t * t)
    rewrite_node(value[0])
    rewrite_node(value[1])
    rewrite_node(value[2])

def rewrite_while(value):
    # While of (t * t)
    rewrite_node(value[0])
    rewrite_node(value[1])

def rewrite_for(value):
    # For of (linno * arg * t * string)
    rewrite_node(value[2])

def rewrite_case(value):
    # Case of (linno * arg * case list)
    (lambda x: x)(value) #no-op

def rewrite_defun(value):
    # Defun of (linno * string * t)
    rewrite_node(value[2])
