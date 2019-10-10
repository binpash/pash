import json
import re

### --- From JSON --- ###

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

def from_standard_json(string):
    string = string.replace("{", "<")
    string = string.replace("}", ">")
    ## TODO: How should I handle tuples?

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

### --- To JSON --- ###

def save_asts_json(asts, json_filename):
    json_string = serialize_asts_to_json(asts)
    with open(json_filename, 'w') as json_file:
        json_file.write(json_string)

def serialize_asts_to_json(asts):
    serialized_asts = [serialize_ast_json(ast) for ast in asts]
    return "\n".join(serialized_asts)

def serialize_ast_json(ast):
    standard_json = json.dumps(ast)
    ocaml_json = from_standard_json(standard_json)
    ## TODO: It seems that it works without reverting all the changes
    ## from standard_json. Is it correct?
    print("Output json:", ocaml_json)
    return ocaml_json
