import json
import re
from ir import *

# json_filename = "../scripts/json/ngrams.sh.json"
json_filename = "../scripts/json/grep.sh.json"

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
        for ast_object in ast_objects:
            # print(json.dumps(ast_object, indent=2))
            print(ast_object)
        return ast_objects

## Checks if the given ASTs are supported
def check_if_asts_supported(ast_objects):
    ## TODO: Implement
    return


### Utils

## TODO: Move to another file
def flatten(l):
    return [item for sublist in l for item in sublist]

## For now these checks are too simple. 
##
## Maybe we can move them to the check_if_ast_is_supported?
def check_pipe(construct, arguments):
    assert(len(arguments) == 2)

def check_command(construct, arguments):
    assert(len(arguments) == 4)

def compile_node(ast_node, nodes):
    # print("Compiling node: {}".format(ast_node))

    construct, arguments = get_kv(ast_node)
    new_nodes = []
    if (construct == 'Pipe'):
        check_pipe(construct, arguments)

        # TODO: Refactor this into a function of its own
        pipe_items = arguments[1]
        new_nodes = flatten([compile_node(pipe_item, []) for pipe_item in pipe_items])
        # Note: I am not sure if it is fine to just give an empty list
        # to each compile node, or whether the compile node needs to
        # get the nodes that were created previously (The concern is
        # whether "compilation" is a map or a fold)
        
        ## TODO: Connect the stdins with the stdouts
        
    elif (construct == 'Command'):
        check_command(construct, arguments)

        command_name = arguments[2][0]
        options = arguments[2][1:]
        new_nodes = [Command(command_name, options=options)]
        
    return new_nodes

## Compiles a given AST to the intermediate representation graph
def compile_ast(ast_object):
    nodes = compile_node(ast_object, [])
    return nodes

## Translation process:
## 1. Parse json to an AST object
## 2. Ensure that the AST is in our supported subset (TODO)
## 3. Compile the AST to out intermediate representation

ast_objects = parse_json_ast(json_filename)
check_if_asts_supported(ast_objects)

for i, ast_object in enumerate(ast_objects):
    print("Compiling AST {}".format(i))
    nodes = compile_ast(ast_object)
    print("Nodes:")
    print(nodes)
