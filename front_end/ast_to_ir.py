import json
import re
from ir import *

# The following files just contain POSIX (no process substitution)
# pipes and commands
#
# json_filename = "../scripts/json/compile.sh.json"
# json_filename = "../scripts/json/grep.sh.json"
# json_filename = "../scripts/json/minimal.sh.json"
# json_filename = "../scripts/json/shortest-scripts.sh.json"
# json_filename = "../scripts/json/spell.sh.json"
# json_filename = "../scripts/json/topn.sh.json"
# json_filename = "../scripts/json/wc.sh.json"
# json_filename = "../scripts/json/wf.sh.json"

# The following contain And operatos together with pipes and commands
#
json_filename = "../scripts/json/ngrams.sh.json" 


# The following is interesting, since it contains command substitution
# (which is parsed as backticks in the Greenberg parser). The
# backticks seem to mean that whatever is in the backticks will
# execute first, and its output will become a string in place of the
# backticks
#
# TODO: Handle appropriately
#
# json_filename = "../scripts/json/page-count.sh.json"

# Unidentified
#
# json_filename = "../scripts/json/maximal.sh.json"


## The json dumper in ocaml seems to print <, >, and parentheses
## instead of {, }, [,]. Therefore we need to replace the characters
## with the correct ones.
def to_standard_json(string):
    string = string.replace("<", "{")
    string = string.replace(">", "}")
    string = string.replace("(", "[")
    string = string.replace(")", "]")

    # After these replacements, single names are written like this:
    # {"Name"} and the parser complains. We just need to remove the
    # braces.
    #
    # Note: I have noticed that the names are always constructors that
    # have no arguments, so they should all be letter characters.
    #
    # Warning: This is not robust at all, but will do for now
    string = re.sub(r'\{\"([A-Za-z]+)\"\}', r'"\1"', string)
    
    return string

## Returns the ast as a object
def parse_json_line(line):
    std_json_line = to_standard_json(line)        
    # print(std_json_line)
    ast_object = json.loads(std_json_line)
    return ast_object

## Returns a list of AST objects
def parse_json_ast(json_filename):
    with open(json_filename) as json_file:
        lines = json_file.readlines()
        ast_objects = [parse_json_line(line) for line in lines]
        # for ast_object in ast_objects:
            # print(json.dumps(ast_object, indent=2))
            # print(ast_object)
        return ast_objects

## Checks if the given ASTs are supported
def check_if_asts_supported(ast_objects):
    ## TODO: Implement
    return


## This checks all ast children nodes of a pipe, and merges the
## consecutive ones to larger IRs.
def combine_pipe(ast_nodes):
    combined_nodes = []
    curr = IR()
    for ast_node in ast_nodes:
        if (isinstance(ast_node, IR)):
            curr.append(ast_node)
        else:
            if (not curr.empty()):
                combined_nodes.append(curr)
                curr = IR()
            
            combined_nodes.append(ast_node)

    if (not curr.empty()):
        combined_nodes.append(curr)
    
    return combined_nodes


## For now these checks are too simple. 
##
## Maybe we can move them to the check_if_ast_is_supported?
def check_pipe(construct, arguments):
    assert(len(arguments) == 2)

def check_command(construct, arguments):
    assert(len(arguments) == 4)

def check_and(construct, arguments):
    assert(len(arguments) == 2)

def compile_arg_char(arg_char):
    key, val = get_kv(arg_char)
    if (key == 'C'):
        return arg_char
    elif (key == 'B'):
        compiled_node = compile_node(val)
        return {key : compiled_node}
    elif (key == 'Q'):
        compiled_val = compile_command_argument(val)
        return {key : compiled_val}
    else:
        ## TODO: Complete this
        return arg_char
    
def compile_command_argument(argument):
    compiled_argument = [compile_arg_char(char) for char in argument]
    return compiled_argument
    
def compile_command_arguments(arguments):
    compiled_arguments = [compile_command_argument(arg) for arg in arguments]
    return compiled_arguments
    
def compile_node(ast_node):
    # print("Compiling node: {}".format(ast_node))

    construct, arguments = get_kv(ast_node)
    
    ## TODO: For simplicity and to get things going, a node that has a
    ## descendant that is backtick should not be distributed for
    ## now.
    ##
    ## Warning: Eventually we should find a way to lift the above assumption
    
    if (construct == 'Pipe'):
        check_pipe(construct, arguments)

        # TODO: Refactor this into a function of its own
        pipe_items = arguments[1]
        compiled_ast = {construct : [arguments[0]] +
                        combine_pipe([compile_node(pipe_item)
                                      for pipe_item in pipe_items])}
        # Note: I am not sure if it is fine to just give an empty list
        # to each compile node, or whether the compile node needs to
        # get the nodes that were created previously (The concern is
        # whether "compilation" is a map or a fold)
        
        ## TODO: Connect the stdins with the stdouts
        
    elif (construct == 'Command'):
        check_command(construct, arguments)

        ## If there are no arguments, the command is just an
        ## assignment
        if(len(arguments[2]) == 0):
            ## TODO: Should we even deal with assignments? If not,
            ## then we can just return them as is
            command_name = []
            options = []
        else:
            command_name = arguments[2][0]
            options = compile_command_arguments(arguments[2][1:])

        ## Question: Should we return the command in an IR if one of
        ## its arguments is a command substitution? Meaning that we
        ## will have to wait for its command to execute first?
        compiled_ast = IR(nodes = [Command(command_name, options=options)])

    elif (construct == 'And'):
        check_and(construct, arguments)
        
        left_node = arguments[0]
        right_node = arguments[1]
        compiled_ast = {construct : [compile_node(left_node), compile_node(right_node)]}
        
    else:
        raise Error("Unimplemented construct: {}".format(construct))
    
    return compiled_ast

## Compiles a given AST to an intermediate representation tree, which
## has some subtrees in it that are graphs representing a distributed
## computation.
##
## The above assumes that subtrees of the AST are disjoint
## computations that can be distributed separately (locally/modularly)
## without knowing about previous or later subtrees that can be
## distributed. Is that reasonable?
def compile_ast(ast_object):
    compiled_ast = compile_node(ast_object)
    return compiled_ast

## Translation process:
## 1. Parse json to an AST object
## 2. Ensure that the AST is in our supported subset (TODO)
## 3. Compile the AST to out intermediate representation

ast_objects = parse_json_ast(json_filename)
check_if_asts_supported(ast_objects)

for i, ast_object in enumerate(ast_objects):
    print("Compiling AST {}".format(i))
    print(ast_object)
    compiled_ast = compile_ast(ast_object)
    print("Compiled AST:")
    print(compiled_ast)
