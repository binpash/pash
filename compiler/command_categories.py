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
    first_opt = format_expanded_arg_chars(options[0])
    if(not (first_opt.startswith("-")
            or first_opt.startswith("s"))
       and ("d" in first_opt
            or "q" in first_opt)):
        return "pure"
    else:
        return "stateless"

def contains_s(option):
    return ((option.startswith("-") and "s" in option)
            or option == "--squeeze-repeats")

def contains_d(option):
    return ((option.startswith("-") and "d" in option)
            or option == "--delete")

def is_tr_pure(options):
    formatted_opts = [format_expanded_arg_chars(option)
                      for option in options]
    set_opts = [opt for opt in formatted_opts if not opt.startswith("-")]
    set1_opt = set_opts[0]
    ## If -s is one of the options (and \n is in the last SET)
    ## If -d is one of the options (and \n is in SET1)
    if((any([contains_s(opt) for opt in formatted_opts])
        and ("\\n" in set_opts[-1]
             or "\\012" in set_opts[-1]))
       or (any([contains_d(opt) for opt in formatted_opts])
           and ("\\n" in set1_opt
                or "\\012" in set1_opt))):
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
    "tr"   : is_tr_pure,
    "uniq" : is_uniq_pure,
}

class Aggregator:
    def __init__(self, aggregator_json):
        self.name = aggregator_json['name']
        self.options = aggregator_json['options']
    
    def __repr__(self):
        return "Aggregator(name={},opts={})".format(self.name, self.options)


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
        # log("inputs-outputs found for:", command_string)
        # log("|--", command_io_from_annotation)
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
        # log("arguments, redirs found for:", command_string)
        # log("|--", command_arguments_redirs)
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

def find_command_properties(command, options):

    command_string = format_arg_chars(command)
    # log("Command to find properties:", command_string)

    assert(isinstance(command_string, str))

    ## Find the properties of the command in the annotation files (if it exists)
    command_properties_from_annotation = get_command_properties_from_annotations(command_string,
                                                                                 options,
                                                                                 config.annotations)
    if (command_properties_from_annotation):
        log("properties:", command_properties_from_annotation, "found for:", command_string)
        return command_properties_from_annotation

    return []

def find_command_aggregator(command, options):

    command_string = format_arg_chars(command)
    # log("Command to find aggregator:", command_string)

    assert(isinstance(command_string, str))

    ## Find the aggregator of the command in the annotation files (if it exists)
    command_aggregator_from_annotation = get_command_aggregator_from_annotations(command_string,
                                                                                 options,
                                                                                 config.annotations)
    if (command_aggregator_from_annotation):
        aggregator_object = Aggregator(command_aggregator_from_annotation)
        log("aggregator:", aggregator_object, "found for:", command_string)
        return aggregator_object

    return None
