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

## The class that contains the aggregator information that was parsed from the annotation
class Aggregator:
    def __init__(self, aggregator_json):
        ## Exactly one of the two should exist in the JSON
        assert('name' in aggregator_json or 'path' in aggregator_json)
        assert(not ('name' in aggregator_json and 'path' in aggregator_json))
        
        ## TODO: Instead of initializing like this, we could keep both, and return the correct information when asked the name.
        if('path' in aggregator_json):
            ## Set the name to be the absolute path
            self.name = os.path.join(config.PASH_TOP, aggregator_json['path'])
        else:
            self.name = aggregator_json['name']
        
        ## By default options are []
        if('options' in aggregator_json):
            self.options = aggregator_json['options']
        else:
            self.options = []
    
    def __repr__(self):
        return "Aggregator(name={},opts={})".format(self.name, self.options)


## The class that contains the mapper information that was mapped from annotation (rare).
class Mapper:
    def __init__(self, mapper_json):
        ## Exactly one of the two should exist in the JSON
        assert('name' in mapper_json or 'path' in mapper_json)
        assert(not ('name' in mapper_json and 'path' in mapper_json))
        
        if('path' in mapper_json):
            ## Set the name to be the absolute path
            self.name = os.path.join(config.PASH_TOP, mapper_json['path'])
        else:
            self.name = mapper_json['name']

        ## By default options are []
        if('options' in mapper_json):
            self.options = mapper_json['options']
        else:
            self.options = []
        
        if ('num_outputs' in mapper_json):
            self.num_outputs = mapper_json['num_outputs']
        else:
            self.num_outputs = 1
    
    def __repr__(self):
        if(self.num_outputs == 1):
            return "Mapper(name={},opts={})".format(self.name, self.options)
        else:
            return "Mapper(name={},opts={},num_outs={})".format(self.name, self.options, self.num_outputs)

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
    log("caruca command to categorize:", command_string)

    assert(isinstance(command_string, str))

    ## Override standard categories
    if (command_string in custom_command_categories):
        log("caruca: Overriding standard category for:", command_string)
        custom_category_fun = custom_command_categories[command_string]
        command_class_from_annotation = custom_category_fun(options)
        # Repeated code is repeated :|
        log("caruca class:", command_class_from_annotation, "found for:", command_string, format_args(options))
        return command_class_from_annotation

    ## Find the class of the command in the annotation files (if it exists)
    command_class_from_annotation = get_command_class_from_annotations(command_string,
                                                                       options,
                                                                       config.annotations)

    log("caruca class:", command_class_from_annotation, "found for:", command_string, format_args(options))

    if (command_class_from_annotation):
        log("class:", command_class_from_annotation, "found for:", command_string)
        return command_class_from_annotation

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

def find_command_mapper_aggregator(command, options):

    command_string = format_arg_chars(command)
    # log("Command to find aggregator:", command_string)

    assert(isinstance(command_string, str))

    ## Find the aggregator of the command in the annotation files (if it exists)
    command_aggregator_from_annotation = get_command_aggregator_from_annotations(command_string,
                                                                                 options,
                                                                                 config.annotations)
    aggregator_object = None
    if (command_aggregator_from_annotation):
        aggregator_object = Aggregator(command_aggregator_from_annotation)
        log("aggregator:", aggregator_object, "found for:", command_string)
        
    ## Find the mapper of the command in the annotation files (if it exists)
    command_mapper_from_annotation = get_command_mapper_from_annotations(command_string,
                                                                         options,
                                                                         config.annotations)
    mapper_object = None
    if (command_mapper_from_annotation):
        mapper_object = Mapper(command_mapper_from_annotation)
        log("mapper:", mapper_object, "found for:", command_string)
    
    return (mapper_object, aggregator_object)

