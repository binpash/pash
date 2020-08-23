from annotations import *
from ir_utils import *

import config

###
### This file contains the logic that determines the category and the
### input outputs of a command. In the future, this should be replaced
### by a parser of the command classification DSL.
###

## TODO: Make the DSL

##
## Custom input output functions for specific commands
##

def cat_input_output(options, stdin, stdout):
    input_stream = [("option", i) for i in range(len(options))]
    return (input_stream, ["stdout"], [])

def comm_input_output(options, stdin, stdout):
    first_opt = format_arg_chars(options[0])
    if(first_opt == "-13"):
        input_opt = format_arg_chars(options[2])
        if(input_opt == "-"):
            in_stream = ["stdin"]
            opt_indices = [("option", i) for i in range(len(options))]
        else:
            in_stream = [("option", 2)]
            opt_indices = [("option", 0), ("option", 1)]
        return (in_stream, ["stdout"], opt_indices)
    elif (first_opt == "-23"):
        input_opt = format_arg_chars(options[1])
        if(input_opt == "-"):
            in_stream = ["stdin"]
            opt_indices = [("option", i) for i in range(len(options))]
        else:
            in_stream = [("option", 1)]
            opt_indices = [("option", 0), ("option", 2)]
        return (in_stream, ["stdout"], opt_indices)
    else:
        assert(false)

## WARNING: This is not complete!! It doesn't handle - for stdin, or
## directories, etc...
##
## TODO: Make this complete
def diff_input_output(options, stdin, stdout):
    in_stream = [("option", i) for i, option in enumerate(options)
                 if not format_arg_chars(option).startswith('-')]
    opt_indices = [("option", i) for i, option in enumerate(options)
                   if format_arg_chars(option).startswith('-')]
    return (in_stream, ["stdout"], opt_indices)

def default_input_output(options, stdin, stdout):
    opt_indices = [("option", i) for i in range(len(options))]
    return (["stdin"], ["stdout"], opt_indices)


##
## Custom category functions for specific commands
##

## TODO: All of these are possibly non-complete
def is_comm_pure(options):
    first_opt = format_arg_chars(options[0])
    if(first_opt == "-13" or first_opt == "-23"):
        return "stateless"
    else:
        return "none"

def is_tr_pure(options):
    first_opt = format_arg_chars(options[0])
    if(first_opt.startswith("-")
       and "d" in first_opt):
        return "pure"
    else:
        return "stateless"

def is_sed_pure(options):
    first_opt = format_arg_chars(options[0])
    if(not first_opt.startswith("-")
       and ("d" in first_opt
            or "q" in first_opt)):
        return "pure"
    else:
        return "stateless"

##
## Dictionaries with the custom functions
##

custom_command_input_outputs = {
    "cat"  : cat_input_output,
    "comm" : comm_input_output,
    "diff" : diff_input_output
}

custom_command_categories = {
    "comm" : is_comm_pure,
    "tr"   : is_tr_pure,
    "sed"  : is_sed_pure,
}


## This function returns the input and output streams of a command.
##
## The input and output lists, contain tuples that refer to options:
## e.g. ("option", 0) or "stdin", "stdout" when they refer to stdin or
## stdout.
##
## At the moment it has just hardcoded knowledge of the inputs and
## outputs of several commands.
##
## By default they are the stdin and the stdout of the node, and they
## are only filled in for commands that we (or the developer) has
## specified a list of input resources that also contains files in the
## arguments.
def find_command_input_output(command, options, stdin, stdout):
    global custom_command_input_outputs

    command_string = format_arg_chars(command)
    # print("Command to categorize:", command_string)

    assert(isinstance(command_string, str))

    if (command_string in custom_command_input_outputs):
        print(" -- Warning: Overriding standard inputs-outputs for:", command_string)
        custom_io_fun = custom_command_input_outputs[command_string]
        return custom_io_fun(options, stdin, stdout)

    ## Find the inputs-outputs of the command in the annotation files (if it exists)
    command_io_from_annotation = get_command_io_from_annotations(command_string,
                                                                 options,
                                                                 config.annotations)
    if (command_io_from_annotation):
        # print("inputs-outputs found for:", command_string)
        return command_io_from_annotation

    return default_input_output(options, stdin, stdout)


## This functions finds and returns a string representing the command category
def find_command_category(command, options):
    global custom_command_categories

    command_string = format_arg_chars(command)
    # print("Command to categorize:", command_string)

    assert(isinstance(command_string, str))

    ## Override standard categories
    if (command_string in custom_command_categories):
        print(" -- Warning: Overriding standard category for:", command_string)
        custom_category_fun = custom_command_categories[command_string]
        return custom_category_fun(options)

    ## Find the class of the command in the annotation files (if it exists)
    command_class_from_annotation = get_command_class_from_annotations(command_string,
                                                                       options,
                                                                       config.annotations)
    if (command_class_from_annotation):
        # print("annotation found for:", command_string)
        return command_class_from_annotation

    # NOTE order of class declaration in definition file is important, as it
    # dictates class precedence in the following search
    for command_class, commands in config.command_classes.items():
        command_list = list(map(get_command_from_definition, commands))

        if (command_string in command_list
            or command_string.split("/")[-1] in command_list):
            return command_class

    return('none')

