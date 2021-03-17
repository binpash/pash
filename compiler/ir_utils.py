### Utils

from util import *

##
## Separate the option from the non-option command arguments
##

def option_args(options):
    option_args_inds = option_args_indices(options)
    args = [option for option, i in option_args_inds]
    return args

def option_args_indices(options):
    non_option_indices = [i for _, i in non_option_args_indices(options)]
    non_option_indices_set = set(non_option_indices)

    ## Find the option args by finding the complement of the
    ## non-option args
    formated_options = [format_arg_chars(opt) for opt in options]
    option_args = [(option, i) for i, option in enumerate(formated_options)
                   if not i in non_option_indices_set]
    return option_args

def non_option_args(options):
    non_option_args_inds = non_option_args_indices(options)
    args = [option for option, i in non_option_args_inds]
    return args

def non_option_args_indices(options):
    formated_options = [format_arg_chars(opt) for opt in options]
    # log(formated_options)

    ## TODO: This might need to become more general
    ##
    ## WARNING: Handling `-` as stdin should not be done for all
    ## commands but only for those that have the stdin-hyphen option.
    args = [(option, i) for i, option in enumerate(formated_options)
            if not option.startswith("-") or option == "-"]
    return args

## This function interleaves option arguments (that might contain Nones)
## with the rest of the arguments
##
## Assumption: rest_arguments does not contain Nones
def interleave_args(opt_arguments, rest_arguments):
    assert(all([arg for arg in rest_arguments if not arg is None]))
    arguments = opt_arguments
    for i in range(len(arguments)):
        if(arguments[i] is None):
            rest_arg = rest_arguments.pop(0)
            arguments[i] = rest_arg
    arguments += rest_arguments
    return arguments

def get_command_from_definition(command_definition):
    if 'command' in command_definition:
        return command_definition['command']

    log('Possible issue with definition file: Missing command in command definition {}'.format(command_definition))
    return ''

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

## These functions check tuple inputs (configuration and streaming ones)
def is_single_input(inputs):
    assert(isinstance(inputs, tuple))
    conf_inputs = inputs[0]
    streaming_inputs = inputs[1]
    return (len(conf_inputs) == 0
            and len(streaming_inputs) == 1)

## This function gets a key and a value from the ast json format
def get_kv(dic):
    return (dic[0], dic[1])

def make_kv(key, val):
    return [key, val]

def string_to_arguments(string):
    return [string_to_argument(word) for word in string.split(" ")]

def string_to_argument(string):
    return [char_to_arg_char(char) for char in string]

## FIXME: This is certainly not complete. It is used to generate the
## AST for the call to the distributed planner. It only handles simple
## characters
def char_to_arg_char(char):
    return ['C' , ord(char)]

def standard_var_ast(string):
    return make_kv("V", ["Normal", False, string, []])

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

def make_command(arguments, redirections=[], assignments = []):
    lineno = 0
    node = make_kv("Command", [lineno, assignments, arguments, redirections])
    return node
