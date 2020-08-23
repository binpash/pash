from ir_utils import *

## Checks if the annotation for that command exists
def get_command_io_from_annotations(command, options, annotations):
    command_ann = get_command_from_annotations(command, options, annotations)
    if(command_ann):
        ## TODO: Interpret those
        inputs = command_ann['inputs']
        outputs = command_ann['outputs']
        return None

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

    ## TODO: Fill in the rest
    return lambda x: False


##
## Helper functions for predicate interpretation
##

def len_args(desired_length, options):
    formated_options = [format_arg_chars(opt) for opt in options]
    # print(formated_options)

    ## filter out the options (starting with `-`)
    args = [opt for opt in formated_options
            if not opt.startswith("-")]
    return (len(args) == desired_length)
