from ir_utils import *

## Checks if the annotation for that command exists
def get_command_io_from_annotations(command, options, annotations):
    command_ann = get_command_from_annotations(command, options, annotations)
    if(command_ann):
        inputs = command_ann['inputs']
        outputs = command_ann['outputs']
        ann_options = []
        if('options' in command_ann):
            ann_options = command_ann['options']
        extracted_inputs = interpret_io_list(inputs, options, ann_options)
        extracted_outputs = interpret_io_list(outputs, options, ann_options)
        option_indices = rest_options(extracted_inputs, extracted_outputs, options)
        return (extracted_inputs, extracted_outputs, option_indices)

def rest_options(inputs, outputs, options):
    input_output_indices = [io[1] for io in inputs + outputs
                            if isinstance(io, tuple)]
    io_indices_set = set(input_output_indices)
    all_indices = [("option", i) for i in range(len(options))
                   if not i in io_indices_set]
    return all_indices

def interpret_io_list(files, options, ann_options):
    io_files = []
    for io in files:
        io_files += interpret_io(io, options, ann_options)
    return io_files

def interpret_io(io, options, ann_options):
    if(io == "stdin"):
        return ["stdin"]
    elif(io == "stdout"):
        return ["stdout"]
    else:
        assert(io.startswith("args"))
        indices = io.split("[")[1].split("]")[0]

        # print(io, options)

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

        ## If the command has the "stdin-hyphen" option turned on,
        ## then it means that `-` should be interpreted as stdin
        if('stdin-hyphen' in ann_options):
            io_list = [handle_stdin_hyphen(io, options) for io in io_list]
        return io_list

def handle_stdin_hyphen(io, options):
    if(isinstance(io, tuple)):
        option = options[io[1]]

        ## TODO: For absolute completeness, not being able to
        ## interpret the argument to `-` means that it should also
        ## bump up the command class (or not be translated to the DFG).
        if(format_arg_chars(option) == '-'):
            return "stdin"
    return io

def get_command_class_from_annotations(command, options, annotations):
    command_ann = get_command_from_annotations(command, options, annotations)
    if(command_ann):
        return command_ann['class']

def get_command_from_annotations(command, options, annotations):
    ## Find an annotation for this command
    relevant_annotations = [ann for ann in annotations if ann['command'] == command]
    if(len(relevant_annotations) == 0):
        return None
    elif(len(relevant_annotations) > 1):
        print("Warning: More than one annotation for command:", command)

    return get_command_from_annotation(command, options, relevant_annotations[0])

def get_command_from_annotation(command, options, annotation):
    assert(annotation['command'] == command)

    cases = annotation['cases']
    case = find_annotation_case(options, cases)
    return case

def find_annotation_case(options, cases):
    for case in cases:
        if(predicate_satisfied(options, case['predicate'])):
            return case

    ## Unreachable
    assert(False)

def predicate_satisfied(options, predicate):
    if(predicate == 'default'):
        return True

    func = interpret_predicate(predicate)
    return func(options)

def interpret_predicate(predicate):
    # print(predicate)
    operator = predicate['operator']
    operands = []
    try:
        operands = predicate['operands']
    except:
        pass
    if(operator == 'len_args_eq'):
        return lambda options: len_args(operands[0], options)
    elif(operator == 'exists'):
        return lambda options: exists_operator(operands, options)
    elif(operator == 'value'):
        return lambda options: value_operator(operands, options)
    elif(operator == 'all'):
        return lambda options: all_operator(operands, options)
    elif(operator == 'or'):
        return lambda options: or_operator(operands, options)
    elif(operator == 'not'):
        return lambda options: not_operator(operands, options)

    ## TODO: Fill in the rest
    return lambda x: print("Uninterpreted operator:", operator); False


##
## Helper functions for predicate interpretation
##

def len_args(desired_length, options):
    args = non_option_args(options)
    return (len(args) == desired_length)

def exists_operator(desired_options, options):
    opt_args_set = set(option_args(options))
    existence = map(lambda opt: opt in opt_args_set, desired_options)
    return any(existence)

## Checks that an option exists and that it has a specific value
def value_operator(option_value, options):
    args_list = format_args(options)
    desired_opt = option_value[0]
    desired_val = option_value[1]
    ## If the desired option does exist, and the next argument is indeed the correct value
    try:
        opt_i = args_list.index(desired_opt)
        val = args_list[opt_i+1]
        return (desired_val == val)
    except:
        return False

def all_operator(desired_options, options):
    opt_args_set = set(option_args(options))
    existence = map(lambda opt: opt in opt_args_set, desired_options)
    return all(existence)

def or_operator(operands, options):
    operand_predicates = map(lambda op: interpret_predicate(op)(options), operands)
    return any(operand_predicates)

## Operands: One predicate
def not_operator(operands, options):
    assert(len(operands) == 1)
    operand = operands[0]
    return not interpret_predicate(operand)(options)

