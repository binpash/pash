from annotations import *
from ir_utils import *
from util import *

import config

###
### This file contains the logic that determines the category and the
### input outputs of a command. In the future, this should be replaced
### by a parser of the command classification DSL.
###

##
## Custom input output functions for specific commands
##


def default_input_output(options):
    opt_indices = [("option", i) for i in range(len(options))]
    return (["stdin"], ["stdout"], opt_indices)


##
## Custom category functions for specific commands
##

## TODO: All of these are possibly non-complete
def is_sed_pure(options):
    first_opt = format_arg_chars(options[0])
    if(not first_opt.startswith("-")
       and ("d" in first_opt
            or "q" in first_opt)):
        return "pure"
    else:
        return "stateless"

## TODO: Move that to annotation
def is_uniq_pure(options):
    if(len(options) > 0):
        first_opt = format_arg_chars(options[0])
        if(first_opt == "-c"):
            return "pure"
    return "parallelizable_pure"


##
## Dictionaries with the custom functions
##

custom_command_input_outputs = {
    
}

custom_command_args_redirs_from_input_outputs = {
    
}

custom_command_categories = {
    "sed"  : is_sed_pure,
    "uniq" : is_uniq_pure,
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
def find_command_input_output(command, options):
    global custom_command_input_outputs

    command_string = format_arg_chars(command)
    # log("Command to categorize:", command_string)

    assert(isinstance(command_string, str))

    if (command_string in custom_command_input_outputs):
        log(" -- Warning: Overriding standard inputs-outputs for:", command_string)
        custom_io_fun = custom_command_input_outputs[command_string]
        return custom_io_fun(options)

    ## Find the inputs-outputs of the command in the annotation files (if it exists)
    command_io_from_annotation = get_command_io_from_annotations(command_string,
                                                                 options,
                                                                 config.annotations)
    if (command_io_from_annotation):
        log("inputs-outputs found for:", command_string)
        log("|--", command_io_from_annotation)
        return command_io_from_annotation

    return default_input_output(options)

## This function is the reverse of the one above. It gives us arguments and redirections
## from inputs and outputs.
def create_command_arguments_redirs(command, options, inputs, outputs):
    global custom_command_args_redirs_from_input_outputs

    command_string = format_arg_chars(command)
    # log("Command to categorize:", command_string)

    assert(isinstance(command_string, str))

    if (command_string in custom_command_args_redirs_from_input_outputs):
        log(" -- Warning: Overriding standard inputs-outputs for:", command_string)
        custom_io_fun = custom_command_args_redirs_from_input_outputs[command_string]
        return custom_io_fun(options)

    ## Find the inputs-outputs of the command in the annotation files (if it exists)
    command_arguments_redirs = construct_args_redirs(command_string,
                                                     options,
                                                     inputs,
                                                     outputs,
                                                     config.annotations)
    if (command_arguments_redirs):
        log("arguments, redirs found for:", command_string)
        log("|--", command_arguments_redirs)
        return command_arguments_redirs

    ## TODO: Implement that
    raise NotImplementedError()
    return default_arguments_redirs(options, inputs, outputs)

## This functions finds and returns a string representing the command category
def find_command_category(command, options):
    global custom_command_categories

    command_string = format_arg_chars(command)
    # log("Command to categorize:", command_string)

    assert(isinstance(command_string, str))

    ## Override standard categories
    if (command_string in custom_command_categories):
        log(" -- Warning: Overriding standard category for:", command_string)
        custom_category_fun = custom_command_categories[command_string]
        return custom_category_fun(options)

    ## Find the class of the command in the annotation files (if it exists)
    command_class_from_annotation = get_command_class_from_annotations(command_string,
                                                                       options,
                                                                       config.annotations)
    if (command_class_from_annotation):
        log("class:", command_class_from_annotation, "found for:", command_string)
        return command_class_from_annotation

    # NOTE order of class declaration in definition file is important, as it
    # dictates class precedence in the following search
    for command_class, commands in config.command_classes.items():
        command_list = list(map(get_command_from_definition, commands))

        if (command_string in command_list
            or command_string.split("/")[-1] in command_list):
            return command_class

    return('none')

