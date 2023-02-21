
from definitions.ast_node import *
from util import *


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



def format_args(args):
    formatted_args = [format_arg_chars(arg_chars) for arg_chars in args]
    return formatted_args

def format_arg_chars(arg_chars):
    chars = [format_arg_char(arg_char) for arg_char in arg_chars]
    return "".join(chars)

##
## BIG TODO: Fix the formating of arg_chars bask to shell scripts and string.
##           We need to do this the proper way using the parser.
##
def format_arg_char(arg_char):
    key, val = get_kv(arg_char)
    if (key == 'C'):
        return str(chr(val))
    elif (key == 'B'):
        # The $() is just for illustration. This is backticks
        return '$({})'.format(val)
    elif (key == 'Q'):
        formated_val = format_arg_chars(val)
        return '"{}"'.format(formated_val)
    elif (key == 'V'):
        return '${{{}}}'.format(val[2])
    elif (key == 'E'):
        ## TODO: This is not right. I think the main reason for the
        ## problems is the differences between bash and the posix
        ## standard.
        # log(" -- escape-debug -- ", val, chr(val))
        non_escape_chars = [92, # \
                            61, # =
                            91, # [
                            93, # ]
                            45, # -
                            58, # :
                            126,# ~
                            42] # *
        if(val in non_escape_chars):
            return '{}'.format(chr(val))
        else:
            return '\{}'.format(chr(val))
    else:
        log("Cannot format arg_char:", arg_char)
        ## TODO: Make this correct
        raise NotImplementedError

## This function finds the first raw character in an argument.
## It needs to be called on an expanded string.
def format_expanded_arg_chars(arg_chars):
    chars = [format_expanded_arg_char(arg_char) for arg_char in arg_chars]
    return "".join(chars)

def format_expanded_arg_char(arg_char):
    key, val = get_kv(arg_char)
    if (key == 'C'):
        return str(chr(val))
    elif (key == 'Q'):
        formated_val = format_expanded_arg_chars(val)
        return '{}'.format(formated_val)
    elif (key == 'E'):
        ## TODO: I am not sure if this should add \ or not
        ##
        ## TODO: This is not right. I think the main reason for the
        ## problems is the differences between bash and the posix
        ## standard.
        # log(" -- escape-debug -- ", val, chr(val))
        non_escape_chars = [92, # \
                            61, # =
                            91, # [
                            93, # ]
                            45, # -
                            58, # :
                            126,# ~
                            42] # *
        if(val in non_escape_chars):
            return '{}'.format(chr(val))
        else:
            return '\{}'.format(chr(val))
    else:
        log("Expanded arg char should not contain:", arg_char)
        ## TODO: Make this correct
        raise ValueError


## TODO: This seems like it should go to ast_util

## This function gets a key and a value from the ast json format
def get_kv(dic):
    return (dic[0], dic[1])

def make_kv(key, val):
    return [key, val]

def string_to_arguments(string):
    return [string_to_argument(word) for word in string.split(" ")]

def string_to_argument(string):
    ret = [char_to_arg_char(char) for char in string]
    return ret

## FIXME: This is certainly not complete. It is used to generate the
## AST for the call to the distributed planner. It only handles simple
## characters
def char_to_arg_char(char):
    return ['C' , ord(char)]

def escaped_char(char):
    return ['E' , ord(char)]

def standard_var_ast(string):
    return make_kv("V", ["Normal", False, string, []])

def make_quoted_variable(string):
    return make_kv("Q", [standard_var_ast(string)])

def quote_arg(arg):
    return make_kv("Q", arg)

def redir_append_stderr_to_string_file(string):
    return make_kv("File",["Append",2,string_to_argument(string)])

def redir_stdout_to_file(arg):
    return make_kv("File",["To", 1, arg])

def redir_file_to_stdin(arg):
    return make_kv("File",["From", 0, arg])

def make_background(body, redirections=[]):
    lineno = 0
    node = make_kv("Background", [lineno, body, redirections])
    return node

def make_backquote(node):
    node = make_kv("B", node)
    return node

def make_subshell(body, redirections=[]):
    lineno = 0
    node = make_kv("Subshell", [lineno, body, redirections])
    return node

def make_command(arguments, redirections=[], assignments=[]):
    lineno = 0
    node = make_kv("Command", [lineno, assignments, arguments, redirections])
    return node

def make_nop():
    return make_command([string_to_argument(":")])

def make_assignment(var, value):
    lineno = 0
    assignment=(var, value)
    assignments=[assignment]
    node = make_kv("Command", [lineno, assignments, [], []])
    return node

def make_semi_sequence(asts):
    if(len(asts) == 0):
        return make_nop()

    if(len(asts) == 1):
        return asts[0]
    else:
        acc = asts[-1]
        ## Remove the last ast
        iter_asts = asts[:-1]
        for ast in iter_asts[::-1]:
            acc = make_kv("Semi", [ast, acc])
        return acc

def make_defun(name, body):
    lineno = 0
    node = make_kv("Defun", [lineno, name, body])
    return node

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