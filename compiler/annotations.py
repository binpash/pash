import config
import os
import json
from ir_utils import *
from util import *

##
## Load annotation files
##

def load_annotation_file(abs_annotation_filename):
    with open(abs_annotation_filename) as annotation_file:
        try:
            annotation = json.load(annotation_file)
            return [annotation]
        except json.JSONDecodeError as err:
            log("WARNING: Could not parse annotation for file:", abs_annotation_filename)
            log("|-- {}".format(err))
            return []

def load_annotation_files(annotation_dir):
    annotations = []
    if(not os.path.isabs(annotation_dir)):
        annotation_dir = os.path.join(config.PASH_TOP, annotation_dir)

    for (dirpath, _dirnames, filenames) in os.walk(annotation_dir):
        json_filenames = [os.path.join(dirpath, filename) for filename in filenames
                          if filename.endswith(".json")]
        curr_annotations = [ann for filename in json_filenames for ann in load_annotation_file(filename) ]
        annotations += curr_annotations
    return annotations
    
## TODO: Make the input/output io_list an object and therefore also define
##       methods that transform inputs/outputs <-> args.

## TODO: Think of a way to create a procedure that:
##       - given an annotation
##       - creates two functions:
##         1. One that takes options and creates inputs, outputs for a DFG Node
##         2. One that takes inputs+outputs of a DFG and creates arguments and redirections
##            for the corresponding command.
##
##       Such a procedure would replace the current ad-hoc functions that 
##       re-read the annotation and then do what they are supposed to do.

##
## Inputs, Outputs -> Command Arguments and redirections
##
def construct_args_redirs(command, options, input_fids, output_fids, annotations):
    command_ann = get_command_from_annotations(command, options, annotations)
    assert(not command_ann is None)
    ## TODO: Move parsing to happen once when the annotation is loaded.
    _, input_args_redirs_assigner_fun = parse_command_inputs_outputs(command_ann['inputs'], True)
    _, output_args_redirs_assigner_fun = parse_command_inputs_outputs(command_ann['outputs'], False)
    ann_options = []
    if('options' in command_ann):
        ann_options = command_ann['options']
    args, redirs = input_args_redirs_assigner_fun(input_fids, ann_options, [], [])
    args, redirs = output_args_redirs_assigner_fun(output_fids, ann_options, args, redirs)
    return args, redirs

def args_redirs_from_io_list(io_list, fids, ann_options, args, redirs):
    rem_fids = fids
    for io in io_list:
        new_args, new_redirs, rem_fids = args_redirs_from_io_list_el(io, rem_fids, ann_options, args, redirs)
        args = new_args
        redirs = new_redirs
    assert(len(rem_fids) == 0)
    return args, redirs

## TODO: We need to handle args[:] followed by stdin by having a look-ahead.
##       It might not be necessary actually. Since it is valid to have them all in
##       args then.
def args_redirs_from_io_list_el(io, fids, ann_options, args, redirs):
    if(io == "stdin"):
        fid = fids[0]
        rem_fids = fids[1:]
        ## Do not redirect if it is from stdin
        if (fid.has_file_descriptor_resource() and fid.resource.is_stdin()):
            return args, redirs, rem_fids
        else:
            new_redirs = redirs
            new_redirs.append(redir_file_to_stdin(fid.to_ast()))
            return args, new_redirs, rem_fids
    elif(io == "stdout"):
        fid = fids[0]
        rem_fids = fids[1:]
        ## Do not redirect if it is for stdout
        if (fid.has_file_descriptor_resource() and fid.resource.is_stdout()):
            return args, redirs, rem_fids
        else:
            new_redirs = redirs
            new_redirs.append(redir_stdout_to_file(fid.to_ast()))
            return args, new_redirs, rem_fids
    else:
        assert(io.startswith("args"))
        indices = io.split("[")[1].split("]")[0]

        if(not ":" in indices):
            index = int(indices)
            ## The argument list is growing with this, so the index might be larger
            args = pad(args, index)
            fid = fids[0]
            rem_fids = fids[1:]
            args[index] = fid
        else:
            start_i_str, end_i_str = indices.split(":")

            start_i = 0
            if(not start_i_str == ""):
                start_i = int(start_i_str)

            ## TODO: We need to handle args[:] followed by stdin by having a look-ahead.

            ## If it has a variable end we need to add all the fids
            end_i = len(fids) + start_i
            if(not end_i_str == ""):
                ## TODO: This might be wrong
                end_i = int(end_i_str)

            ## The argument list is growing with this, so the index might be larger
            args = pad(args, end_i - 1)

            ## All the arguments in the required range must be None
            assert(len([arg for arg in args[start_i:end_i] if not arg is None]) == 0)
            for i in range(start_i,end_i):
                args[i] = fids[i-start_i]

            ## Remove the used fids
            rem_fids = fids[(end_i-start_i):]

        ## If the command has the "stdin-hyphen" option turned on,
        ## then it means that `-` should be interpreted as stdin
        ## TODO: What should we do for stdin-hyphen?
        # if('stdin-hyphen' in ann_options):
        #     io_list = [handle_stdin_hyphen(io, options) for io in io_list]
        return args, redirs, rem_fids

## This function parses the command inputs and creates a function 
## that can be extracts inputs from options
##
## TODO: Come up with the conditions that need to hold for these functions
##
## TODO: Unify both of these creations into one. Now there is a lot of duplication.
def parse_command_inputs_outputs(inputs_outputs, parsing_inputs : bool):
    if(isinstance(inputs_outputs, list)):
        input_assigner_fun = lambda options, ann_options: interpret_io_list(inputs_outputs, options, ann_options, parsing_inputs)
        args_redirs_assigner_fun = lambda fids, ann_options, args, redirs: args_redirs_from_io_list(inputs_outputs,
                                                                                                    fids,
                                                                                                    ann_options,
                                                                                                    args,
                                                                                                    redirs)
    else:
        configuration_inputs = inputs_outputs["configuration"]
        standard_inputs = inputs_outputs["standard"]
        input_assigner_fun = lambda options, ann_options: assign_configuration_standard_inputs(configuration_inputs, standard_inputs, options, ann_options, parsing_inputs)
        all_inputs = configuration_inputs + standard_inputs
        args_redirs_assigner_fun = lambda fids, ann_options, args, redirs: args_redirs_from_io_list(all_inputs,
                                                                                                    fids,
                                                                                                    ann_options,
                                                                                                    args,
                                                                                                    redirs)

    return (input_assigner_fun, args_redirs_assigner_fun)
    
def assign_configuration_standard_inputs(configuration_inputs, standard_inputs, options, ann_options, parsing_inputs):
    extracted_config_inputs, options_to_rem1 = interpret_io_list(configuration_inputs, options, ann_options, parsing_inputs)
    extracted_standard_inputs, options_to_rem2 = interpret_io_list(standard_inputs, options, ann_options, parsing_inputs)
    extracted_inputs = (extracted_config_inputs, extracted_standard_inputs)
    return (extracted_inputs, options_to_rem1 + options_to_rem2)


## Checks if the annotation for that command exists
def get_command_io_from_annotations(command, options, annotations):
    command_ann = get_command_from_annotations(command, options, annotations)
    if(command_ann):
        log("Annotation:", command_ann)
        ## TODO: Move parsing to happen once when the annotation is loaded.
        input_assigner_fun, _ = parse_command_inputs_outputs(command_ann['inputs'], True)
        output_assigner_fun, _ = parse_command_inputs_outputs(command_ann['outputs'], False)
        ann_options = []
        if('options' in command_ann):
            ann_options = command_ann['options']
        extracted_inputs, options_to_rem1 = input_assigner_fun(options, ann_options)
        extracted_outputs, options_to_rem2 = output_assigner_fun(options, ann_options)

        ## Some options do not have to be considered for option_indices
        ## At the moment this is only `-` for stdin-hyphen
        if(isinstance(extracted_inputs, tuple)):
            options_to_ignore = options_to_rem1 + options_to_rem2 + extracted_inputs[0] + extracted_inputs[1] + extracted_outputs
        else:
            options_to_ignore = options_to_rem1 + options_to_rem2 + extracted_inputs + extracted_outputs
        option_indices = rest_options(options_to_ignore, options)
        return (extracted_inputs, extracted_outputs, option_indices)

def rest_options(options_to_ignore, options):
    input_output_indices = [io[1] for io in options_to_ignore
                            if isinstance(io, tuple)]
    io_indices_set = set(input_output_indices)
    all_indices = [("option", i) for i in range(len(options))
                   if not i in io_indices_set]
    return all_indices

def interpret_io_list(files, options, ann_options, parsing_inputs):
    log("Files:", files)
    io_files = []
    options_to_remove = []
    for io in files:
        new_io_files, new_options_to_remove = interpret_io(io, options, ann_options)
        io_files += new_io_files
        options_to_remove += new_options_to_remove
    ## If there is an annotation option "empty-args-stdin" 
    ##   and there is no io_file until now, we need to put stdin
    if(len(new_io_files) == 0
       and parsing_inputs
       and 'empty-args-stdin' in ann_options):
        assert(len(options_to_remove) == 0)
        io_files = ["stdin"]
    return (io_files, options_to_remove)

def interpret_io(io, options, ann_options):
    if(io == "stdin"):
        return (["stdin"], [])
    elif(io == "stdout"):
        return (["stdout"], [])
    else:
        assert(io.startswith("args"))
        indices = io.split("[")[1].split("]")[0]

        # log(io, options)

        ## Find the file arguments (and their indices)
        args_indices = non_option_args_indices(options)

        ## Single index
        if(not ":" in indices):
            index = int(indices)
            io_list = [("option", args_indices[index][1])]
        else:
            start_i_str, end_i_str = indices.split(":")

            start_i = 0
            if(not start_i_str == ""):
                start_i = int(start_i_str)

            end_i = len(args_indices)
            if(not end_i_str == ""):
                end_i = int(end_i_str)

            io_list = []
            for _, i in args_indices[start_i:end_i]:
                io_list.append(("option", i))

        options_to_remove = []
        ## If the command has the "stdin-hyphen" option turned on,
        ## then it means that `-` should be interpreted as stdin
        if('stdin-hyphen' in ann_options):
            ## We need to remove the `-` from the options.
            ## TODO: This is not a complete solution. Normally we should
            ##       keep this info in the DFGNode somewhere. The compilation
            ##       Command <-> DFGNode is not a bijection at the moment.
            io_list, options_to_remove = handle_stdin_hyphen(io_list, options)
        return (io_list, options_to_remove)

def handle_stdin_hyphen(io_list, options):
    options_to_remove = []
    for i in range(len(io_list)):
        io = io_list[i]
        if(isinstance(io, tuple)):
            option = options[io[1]]

            ## TODO: For absolute completeness, not being able to
            ## interpret the argument to `-` means that it should also
            ## bump up the command class (or not be translated to the DFG).
            if(format_arg_chars(option) == '-'):
                io_list[i] = "stdin"
                ## Remove `-` from options
                options_to_remove.append(io)
    return (io_list, options_to_remove)


def get_command_class_from_annotations(command, options, annotations):
    command_ann = get_command_from_annotations(command, options, annotations)
    if(command_ann):
        return command_ann['class']

## TODO: Find a general way to handle arbitrary paths etc
def get_command_from_annotations(command_path, options, annotations):
    ## Get only the basename from the path of the command
    command_basename = os.path.basename(command_path)
    ## Find an annotation for this command
    relevant_annotations = [ann for ann in annotations if ann['command'] == command_basename]
    if(len(relevant_annotations) == 0):
        return None
    elif(len(relevant_annotations) > 1):
        log("Warning: More than one annotation for command:", command_basename)

    return get_command_from_annotation(command_basename, options, relevant_annotations[0])

def get_command_from_annotation(command, options, annotation):
    assert(annotation['command'] == command)

    cases = annotation['cases']
    ## Maps short to long options
    short_long_correspondence = get_short_long_mapping(annotation)
    case = find_annotation_case(options, short_long_correspondence, cases)

    ## Extend the case with the options from the main annotation
    if('options' in annotation):
        general_options = annotation['options']
        case_options = []
        if('options' in case):
            case_options = case['options']
        case['options'] = general_options + case_options

    return case

## TODO: Get the short-long mapping and pass it to predicates
def get_short_long_mapping(annotation):
    short_long_correspondence = {}
    if('short-long' in annotation):
        short_long_json = annotation['short-long']
        # log(short_long_json)
        for pair in short_long_json:
            short_opt=pair['short']
            long_opt=pair['long']
            short_long_correspondence[short_opt] = long_opt

    return short_long_correspondence

## Uses the short-long correspondence to map as many short options
## to long.
def transform_opts_to_long(opts, short_long):
    # log("Options to transform:", opts)
    new_opts = opts[:]
    for i in range(len(new_opts)):
        opt = new_opts[i]
        if(opt in short_long):
            new_opts[i] = short_long[opt]
    return new_opts

def find_annotation_case(options, short_long, cases):
    for case in cases:
        if(predicate_satisfied(options, short_long, case['predicate'])):
            return case

    ## Unreachable
    assert(False)

def predicate_satisfied(options, short_long, predicate):
    if(predicate == 'default'):
        return True

    func = interpret_predicate(predicate, short_long)
    return func(options)

## TODO: Handle val_opt_eq

def interpret_predicate(predicate, short_long):
    # log(predicate)
    operator = predicate['operator']
    operands = []
    try:
        operands = predicate['operands']
    except:
        pass
    if(operator == 'len_args_eq'):
        return lambda options: len_args(operands[0], options)
    elif(operator == 'exists'):
        return lambda options: exists_operator(operands, short_long, options)
    elif(operator == 'value'):
        return lambda options: value_operator(operands, short_long, options)
    elif(operator == 'all'):
        return lambda options: all_operator(operands, short_long, options)
    elif(operator == 'or'):
        return lambda options: or_operator(operands, short_long, options)
    elif(operator == 'not'):
        return lambda options: not_operator(operands, short_long, options)

    ## TODO: Fill in the rest
    return lambda x: log("Uninterpreted operator:", operator); False


##
## Helper functions for predicate interpretation
##

def len_args(desired_length, options):
    args = non_option_args(options)
    return (len(args) == desired_length)

def exists_operator(desired_options, short_long, options):
    ## First get all option args (ones that start with '-')
    all_opt_args = option_args(options)
    ## Then transform the to their long variants (if they exist)
    long_opt_args = transform_opts_to_long(all_opt_args, short_long)
    opt_args_set = set(long_opt_args)
    ## Also transform the desired options to their long variants (if they exist)
    long_desired_options = transform_opts_to_long(desired_options, short_long)
    existence = map(lambda opt: opt in opt_args_set, long_desired_options)
    return any(existence)

## Checks that an option exists and that it has a specific value
def value_operator(option_value, short_long, options):
    ## First, translate the options to strings
    args_list = format_args(options)
    ## Then transform any options to their long variants (if they exist)
    ## TODO: Warning, this could also replace non-option arguments with long
    long_args_list = transform_opts_to_long(args_list, short_long)
    desired_opt = option_value[0]
    desired_opt_long = transform_opts_to_long([desired_opt], short_long)[0]
    desired_val = option_value[1]
    ## If the desired option does exist, and the next argument is indeed the correct value
    try:
        opt_i = long_args_list.index(desired_opt_long)
        val = args_list[opt_i+1]
        return (desired_val == val)
    except:
        return False

## This is a shortcut for and . exists. It means that all
## desired options should exist in options.
def all_operator(desired_options, short_long, options):
    ## First get all option args (ones that start with '-')
    all_opt_args = option_args(options)
    ## Then transform the to their long variants (if they exist)
    long_opt_args = transform_opts_to_long(all_opt_args, short_long)
    opt_args_set = set(long_opt_args)
    ## Also transform the desired options to their long variants (if they exist)
    long_desired_options = transform_opts_to_long(desired_options, short_long)
    existence = map(lambda opt: opt in opt_args_set, long_desired_options)
    return all(existence)

def or_operator(operands, short_long, options):
    operand_predicates = map(lambda op: interpret_predicate(op, short_long)(options), operands)
    return any(operand_predicates)

## Operands: One predicate
def not_operator(operands, short_long, options):
    assert(len(operands) == 1)
    operand = operands[0]
    return not interpret_predicate(operand, short_long)(options)

