from definitions.ast_node import *
## TODO: These calls need to be moved here, they don't make sense in ir_utils
from ir_utils import string_to_argument, make_kv

## This class is used by the preprocessor in ast_to_ir
class PreprocessedAST:
    def __init__(self, ast, replace_whole, non_maximal, something_replaced=True, last_ast=False):
        assert(isinstance(ast, AstNode))
        self.ast = ast
        self.replace_whole = replace_whole
        self.non_maximal = non_maximal
        self.something_replaced = something_replaced
        self.last_ast = last_ast

    def should_replace_whole_ast(self):
        return self.replace_whole

    def is_non_maximal(self):
        return self.non_maximal
    
    def will_anything_be_replaced(self):
        return self.something_replaced

    def is_last_ast(self):
        return self.last_ast

## This class represents text that was not modified at all by preprocessing, and therefore does not
## need to be unparsed.
class UnparsedScript:
    def __init__(self, text):
        self.text = text


##
## Pattern matching for the AST
##

def check_if_ast_is_supported(construct, arguments, **kwargs):
    return

def ast_match_untyped(untyped_ast_object, cases, *args):
    ## TODO: This should construct the complete AstNode object (not just the surface level)
    ## TODO: Remove this and then at some point make real proper use of the AstNode
    ast_node = AstNode(untyped_ast_object)
    if ast_node.construct is AstNodeConstructor.PIPE:
        ast_node.check(children_count = lambda : len(ast_node.items) >= 2)
    return ast_match(ast_node, cases, *args)

def ast_match(ast_node, cases, *args):
    ## TODO: Remove that once `ast_match_untyped` is fixed to
    ##       construct the whole AstNode object.
    if(not isinstance(ast_node, AstNode)):
        return ast_match_untyped(ast_node, cases, *args)

    return cases[ast_node.construct.value](*args)(ast_node)

##
## Make some nodes
##

def make_echo_ast(argument, var_file_path):
    nodes = []
    ## Source variables if present
    if(not var_file_path is None):
        arguments = [string_to_argument("source"), string_to_argument(var_file_path)]

        line_number = 0
        node = make_kv('Command', [line_number, [], arguments, []])
        nodes.append(node)

    ## Reset the exit status
    variable_arg = make_kv('V', ['Normal', "false", 'pash_previous_exit_status', []])
    arguments = [string_to_argument("exit"), [variable_arg]]
    exit_node = make_kv('Command', [0, [], arguments, []])
    node = make_kv('Subshell', [0, exit_node, []])
    nodes.append(node)

    ## Reset the input arguments
    variable_arg = make_kv('V', ['Normal', "false", 'pash_input_args', []])
    arguments = [string_to_argument("set"), string_to_argument("--"), [variable_arg]]
    set_node = make_kv('Command', [0, [], arguments, []])
    nodes.append(set_node)

    arguments = [string_to_argument("echo"), string_to_argument("-n"), argument]

    line_number = 0
    node = make_kv('Command', [line_number, [], arguments, []])
    nodes.append(node)
    return nodes