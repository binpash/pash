import json
import re

json_filename = "../scripts/json/ngrams.sh.json"

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
            print(json.dumps(ast_object, indent=2))
            print(ast_object)
        return ast_objects

## Checks if the given ASTs are supported
def check_if_asts_supported(ast_objects):
    ## TODO: Implement
    return
    
## Translation process:
## 1. Parse json to an AST object
## 2. Ensure that the AST is in our supported subset (TODO)
## 3. Compile the AST to out intermediate representation

ast_objects = parse_json_ast(json_filename)

check_if_asts_supported(ast_objects)
