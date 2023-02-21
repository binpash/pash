import json
import config
from shell_ast.ast_node import CustomJSONEncoder
from subprocess import run, PIPE

from util import *

### --- From JSON --- ###

## Returns the ast as a object
def parse_json_line(json_line):
    ast_object = json.loads(json_line)
    return ast_object

def parse_json_ast_string(json_string):
    ## Solve a bug where an empty string leads to a json decoder error
    if(json_string == ""):
        log("Warning: Empty json string was given for decoding.")
        return []
    stripped_json_string = json_string.strip()
    lines = stripped_json_string.split("\n")
    ast_objects = [parse_json_line(line) for line in lines]
    return ast_objects

## Returns a list of AST objects
def parse_json_ast(json_filename):
    with open(json_filename) as json_file:
        file_string = json_file.read()
        return parse_json_ast_string(file_string)

### --- To JSON --- ###

def save_asts_json(asts, json_filename):
    json_string = serialize_asts_to_json(asts)
    with open(json_filename, 'w') as json_file:
        json_file.write(json_string)

def serialize_asts_to_json(asts):
    serialized_asts = [serialize_ast_json(ast) for ast in asts]
    return "\n".join(serialized_asts)

def serialize_ast_json(ast):
    standard_json = json.dumps(ast, cls=CustomJSONEncoder)
    return standard_json

### --- AST to Shell --- ###


def json_to_shell(json_string):
    subproc = run([config.PRINTER_BINARY], stdout=PIPE, input=json_string,
                  encoding='ascii', check=True)
    return subproc.stdout

def ast_to_shell(ast, verbose=True):
    ast_json = serialize_ast_json(ast)
    if verbose:
        print(ast_json)
    shell_string = json_to_shell(ast_json)
    return shell_string

