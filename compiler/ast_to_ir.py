from ir import *
from union_find import *
from definitions.ast_node import *
from definitions.ast_node_c import *
from util import *
from json_ast import save_asts_json
from parse import parse_shell, from_ir_to_shell
import subprocess

import config

import pickle


##
## Compile AST -> Extended AST with IRs
##

## Compiles a given AST to an intermediate representation tree, which
## has some subtrees in it that are graphs representing a distributed
## computation.
##
## The above assumes that subtrees of the AST are disjoint
## computations that can be distributed separately (locally/modularly)
## without knowing about previous or later subtrees that can be
## distributed. Is that reasonable?
compile_cases = {
        "Pipe": (lambda fileIdGen, config:
                 lambda ast_node: compile_node_pipe(ast_node, fileIdGen, config)),
        "Command": (lambda fileIdGen, config:
                    lambda ast_node: compile_node_command(ast_node, fileIdGen, config)),
        "And": (lambda fileIdGen, config:
                lambda ast_node: compile_node_and_or_semi(ast_node, fileIdGen, config)),
        "Or": (lambda fileIdGen, config:
               lambda ast_node: compile_node_and_or_semi(ast_node, fileIdGen, config)),
        "Semi": (lambda fileIdGen, config:
                 lambda ast_node: compile_node_and_or_semi(ast_node, fileIdGen, config)),
        "Redir": (lambda fileIdGen, config:
                  lambda ast_node: compile_node_redir_subshell(ast_node, fileIdGen, config)),
        "Subshell": (lambda fileIdGen, config:
                     lambda ast_node: compile_node_redir_subshell(ast_node, fileIdGen, config)),
        "Background": (lambda fileIdGen, config:
                       lambda ast_node: compile_node_background(ast_node, fileIdGen, config)),
        "Defun": (lambda fileIdGen, config:
                  lambda ast_node: compile_node_defun(ast_node, fileIdGen, config)),
        "For": (lambda fileIdGen, config:
                  lambda ast_node: compile_node_for(ast_node, fileIdGen, config))
        }

preprocess_cases = {
    "Pipe": (lambda irFileGen, config:
             lambda ast_node: preprocess_node_pipe(ast_node, irFileGen, config)),
    "Command": (lambda irFileGen, config:
                lambda ast_node: preprocess_node_command(ast_node, irFileGen, config)),
    "For": (lambda irFileGen, config:
            lambda ast_node: preprocess_node_for(ast_node, irFileGen, config))
}

ir_cases = {
        ## Note: We should never encounter a Pipe construct, since all
        ## of them must have been become IRs
        ##
        ## "Pipe": (lambda c, b, i: {c : [b, i]}),
        "Command":   (lambda irFileGen, config:
                      lambda ast_node: replace_irs_command(ast_node, irFileGen, config)),
        "And" :      (lambda irFileGen, config:
                      lambda ast_node: replace_irs_and_or_semi(ast_node, irFileGen, config)),
        "Or" :       (lambda irFileGen, config:
                      lambda ast_node: replace_irs_and_or_semi(ast_node, irFileGen, config)),
        "Semi" :     (lambda irFileGen, config:
                      lambda ast_node: replace_irs_and_or_semi(ast_node, irFileGen, config)),

        ## TODO: Complete these
        # "Redir" : (lambda c, l, n, r: compile_node_redir(c, l, n, r, fileIdGen)),
        # "Subshell" : (lambda c, l, n, r: compile_node_subshell(c, l, n, r, fileIdGen)),
        # "Background" : (lambda c, l, n, r: compile_node_background(c, l, n, r, fileIdGen)),
        # "Defun" : (lambda c, l, n, b: compile_node_defun(c, l, n, b, fileIdGen)),
        "For" :      (lambda irFileGen, config:
                      lambda ast_node: replace_irs_for(ast_node, irFileGen, config)),
        }

def compile_asts(ast_objects, fileIdGen, config):
    compiled_asts = []
    acc_ir = None
    for i, ast_object in enumerate(ast_objects):
        # print("Compiling AST {}".format(i))
        # print(ast_object)

        ## Compile subtrees of the AST to out intermediate representation
        compiled_ast = compile_node(ast_object, fileIdGen, config)

        # print("Compiled AST:")
        # print(compiled_ast)

        ## If the accumulator contains an IR (meaning that the
        ## previous commands where run in background), union it with
        ## the current returned ast.
        if (not acc_ir is None):

            if (isinstance(compiled_ast, IR)):
                acc_ir.union(compiled_ast)
            else:
                ## TODO: Make this union the compiled_ast with the
                ## accumulated IR, since the user wanted to run these
                ## commands in parallel (Is that correct?)
                # acc_ir.union(IR([compiled_ast]))
                compiled_asts.append(acc_ir)
                acc_it = None
                compiled_asts.append(compiled_ast)

            ## If the current compiled ast not in background (and so
            ## the union isn't in background too), stop accumulating
            if (not acc_ir is None
                and not acc_ir.is_in_background()):
                compiled_asts.append(acc_ir)
                acc_ir = None
        else:
            ## If the compiled ast is in background, start
            ## accumulating it
            if (isinstance(compiled_ast, IR)
                and compiled_ast.is_in_background()):
                acc_ir = compiled_ast
            else:
                compiled_asts.append(compiled_ast)

    ## The final accumulator
    if (not acc_ir is None):
        compiled_asts.append(acc_ir)

    return compiled_asts


def compile_node(ast_object, fileIdGen, config):
    global compile_cases
    return ast_match(ast_object, compile_cases, fileIdGen, config)

def compile_node_pipe(ast_node, fileIdGen, config):
    compiled_pipe_nodes = combine_pipe([compile_node(pipe_item, fileIdGen, config)
                                        for pipe_item in ast_node.items])

    if (len(compiled_pipe_nodes) == 1):
        ## Note: When calling combine_pipe_nodes (which
        ##       optimistically distributes all the children of a
        ##       pipeline) the compiled_pipe_nodes should always
        ##       be one IR
        compiled_ir = compiled_pipe_nodes[0]
        ## Note: Save the old ast for the end-to-end prototype
        old_ast_node = make_kv(ast_node.construct.value, [ast_node.is_background, ast_node.items])
        compiled_ir.set_ast(old_ast_node)
        ## Set the IR background so that it can be parallelized with
        ## the next command if the pipeline was run in background
        compiled_ir.set_background(ast_node.is_background)
        ## TODO: If the pipeline is in background, I also have to
        ## redirect its stdin, stdout
        compiled_ast = compiled_ir
    else:
        ## Note: This should be unreachable with the current combine pipe
        assert(False)
        compiled_ast = make_kv(construct, [arguments[0]] + [compiled_pipe_nodes])
    return compiled_ast

## This combines all the children of the Pipeline to an IR, even
## though they might not be IRs themselves. This means that an IR
## might contain mixed commands and ASTs. The ASTs can be
## (conservatively) considered as stateful commands by default).
def combine_pipe(ast_nodes):
    ## Initialize the IR with the first node in the Pipe
    if (isinstance(ast_nodes[0], IR)):
        combined_nodes = ast_nodes[0]
    else:
        ## FIXME: This one will not work. The IR of an AST node
        ##        doesn't have any stdin or stdout.
        combined_nodes = IR([ast_nodes[0]])

    ## Combine the rest of the nodes
    for ast_node in ast_nodes[1:]:
        if (isinstance(ast_node, IR)):
            combined_nodes.pipe_append(ast_node)
        else:
            ## FIXME: This one will not work. The IR of an AST node
            ##        doesn't have any stdin or stdout.
            combined_nodes.pipe_append(IR([ast_node]))

    return [combined_nodes]

def compile_node_command(ast_node, fileIdGen, config):
    construct_str = ast_node.construct.value
    old_ast_node = make_kv(construct_str, [ast_node.line_number,
        ast_node.assignments, ast_node.arguments, ast_node.redir_list])

    ## TODO: Do we need the line number?

    ## Compile assignments and redirection list
    compiled_assignments = compile_assignments(ast_node.assignments, fileIdGen, config)
    compiled_redirections = compile_redirections(ast_node.redir_list, fileIdGen, config)

    ## If there are no arguments, the command is just an
    ## assignment
    if(len(ast_node.arguments) == 0):
        ## Just compile the assignments. Specifically compile the
        ## assigned values, because they might have command
        ## substitutions etc..
        compiled_ast = make_kv(construct_str, [ast_node.line_number] +
                               [compiled_assignments] + [ast_node.arguments, compiled_redirections])
    else:
        arguments = ast_node.arguments
        command_name = arguments[0]
        options = compile_command_arguments(arguments[1:], fileIdGen, config)

        stdin_fid = fileIdGen.next_file_id()
        stdout_fid = fileIdGen.next_file_id()
        ## Question: Should we return the command in an IR if one of
        ## its arguments is a command substitution? Meaning that we
        ## will have to wait for its command to execute first?
        ##
        ## ANSWER: Kind of. If a command has a command substitution or
        ## anything that evaluates we should add it to the IR, but we
        ## should also make sure that its category is set to the most
        ## general one. That means that it can be executed
        ## concurrently with other commands, but it cannot be
        ## parallelized.
        command = create_command_assign_file_identifiers(old_ast_node, fileIdGen,
                                                         command_name, options,
                                                         stdin=stdin_fid, stdout=stdout_fid,
                                                         redirections=compiled_redirections)

        ## Don't put the command in an IR if it is creates some effect
        ## (not stateless or pure)
        if (command.is_at_most_pure()):
            compiled_ast = IR([command],
                              stdin = [stdin_fid],
                              stdout = [stdout_fid])
            compiled_ast.set_ast(old_ast_node)
        else:
            compiled_arguments = compile_command_arguments(arguments, fileIdGen, config)
            compiled_ast = make_kv(construct_str,
                                   [ast_node.line_number, compiled_assignments,
                                    compiled_arguments, compiled_redirections])

    return compiled_ast

def compile_node_and_or_semi(ast_node, fileIdGen, config):
    compiled_ast = make_kv(ast_node.construct.value,
            [compile_node(ast_node.left_operand, fileIdGen, config),
             compile_node(ast_node.right_operand, fileIdGen, config)])
    return compiled_ast

def compile_node_redir_subshell(ast_node, fileIdGen, config):
    compiled_node = compile_node(ast_node.node, fileIdGen, config)

    if (isinstance(compiled_node, IR)):
        ## TODO: I should use the redir list to redirect the files of
        ##       the IR accordingly
        compiled_ast = compiled_node
    else:
        compiled_ast = make_kv(ast_node.construct.value, [ast_node.line_number,
            compiled_node, ast_node.redir_list])

    return compiled_ast

def compile_node_background(ast_node, fileIdGen, config):
    compiled_node = compile_node(ast_node.node, fileIdGen, config)

    ## TODO: I should use the redir list to redirect the files of
    ##       the IR accordingly
    if (isinstance(compiled_node, IR)):
        ## TODO: Redirect the stdout, stdin accordingly
        compiled_node.set_background(True)
        compiled_ast = compiled_node
    else:
        ## Note: It seems that background nodes can be added in
        ##       the distributed graph similarly to the children
        ##       of pipelines.
        ##
        ## Question: What happens with the stdin, stdout. Should
        ## they be closed?
        compiled_ast = IR([compiled_node])

    return compiled_ast

def compile_node_defun(ast_node, fileIdGen, config):
    ## It is not clear how we should handle functions.
    ##
    ## - Should we transform their body to IR?
    ## - Should we handle calls to the functions as commands?
    ##
    ## It seems that we should do both. But we have to think if
    ## this introduces any possible problem.

    ## TODO: Investigate whether it is fine to just compile the
    ##       body of functions.
    compiled_body = compile_node(ast_node.body, fileIdGen, config)
    return make_kv(construct, [ast_node.line_number, ast_node.name, compiled_body])

def compile_node_for(ast_node, fileIdGen, config):
    ## TODO: Investigate what kind of check could we do to make a for
    ## loop parallel
    compiled_ast = make_kv(ast_node.construct.value,
                           [ast_node.line_number,
                            compile_command_argument(ast_node.argument, fileIdGen, config),
                            compile_node(ast_node.body, fileIdGen, config),
                            ast_node.variable])
    return compiled_ast


## This function checks if a word is safe to expand (i.e. if it will 
## not have unpleasant side-effects)
def safe_to_expand(arg_char):
    key, val = get_kv(arg_char)
    if (key in ['V']): # Variable
        return True
    return False

def safe_arg(arg):
    return all([safe_arg_char(arg_char) for arg_char in arg])

def safe_arg_char(arg_char):
    key, val = get_kv(arg_char)
    # character, escaped---noop, but safe
    if (key in ['C', 'E']): 
        return True
    # tilde --- only reads system state, safe to do early assuming no writes to HOME prior
    elif (key in ['T']):
        return True # TODO 2020-11-24 MMG modified variable set? take in/output written vars...
    # arithmetic -- depends on what we have
    elif (key == 'A'):
        return safe_arith(val)
    # quoted -- safe if its contents are safe
    elif (key == 'Q'):
        return safe_arg(key)
    # variables -- safe if the format is safe as are the remaining words
    elif (key == 'V'):
        return True # TODO depends on format...
    # command substitution -- depends on the command
    elif (key == 'B'):
        return safe_command(val)
    
    raise ValueError("bad key {}, expected one of CETAVQB".format(key))

def safe_var(fmt, null, var, arg):
    if (fmt in ['Normal', 'Length']):
        return True
    elif (fmt in ['Minus', 'Plus', 'Question', 'TrimR', 'TrimRMax', 'TrimL', 'TrimLMax']):
        return safe_arg(arg)
    elif (fmt in ['Assign']):
        return False # TODO 2020-11-24 MMG unless we know `var` is set

    raise ValueError("bad parameter format {}".format(fmt))

def safe_arith(arg):
    # operations are safe
    # `+=` and `=` and family are UNSAFE
    # NONPOSIX: `++` and `--` are UNSAFE
    # `op="+=1"; $((x $op))` is UNSAFE

    # TODO 2020-11-24 MMG
    # to determine safety, we:
    #   (a) check that every arg_char here is safe
    #   (b) pre-parse it symbolically well enough to ensure that no mutating operations occur
    return False

def safe_command(command):
    # TODO 2020-11-24 MMG which commands are safe to run in advance?
    # TODO 2020-11-24 MMG how do we differentiate it being safe to do nested expansions?
    return False

def make_echo_ast(arg_char):
    arguments = [string_to_argument("echo"), string_to_argument("-n"), [arg_char]]

    line_number = 0
    node = make_kv('Command', [line_number, [], arguments, []])
    return node

## TODO: Move this function somewhere more general
def execute_shell_asts(asts):
    ir_filename = os.path.join("/tmp", get_random_string())
    save_asts_json(asts, ir_filename)
    output_script = from_ir_to_shell(ir_filename)
    # print(output_script)
    exec_obj = subprocess.run(["/bin/bash"], input=output_script, 
                              capture_output=True,
                              text=True)
    exec_obj.check_returncode()
    # print(exec_obj.stdout)
    return exec_obj.stdout

## TODO: Properly parse the output of the shell script
def parse_string_to_arg_char(arg_char_string):
    # print(arg_char_string)
    return ['Q', string_to_argument(arg_char_string)]

def naive_expand(arg_char, config):
    ## Create an AST node that "echo"s the argument
    echo_ast = make_echo_ast(arg_char)

    ## Execute the echo AST by unparsing it to shell
    ## and calling bash
    expanded_string = execute_shell_asts([echo_ast])

    ## Parse the expanded string back to an arg_char
    expanded_arg_char = parse_string_to_arg_char(expanded_string)
    
    ## TODO: Handle any errors
    # print(expanded_arg_char)
    return expanded_arg_char



## This function expands an arg_char. 
## At the moment it is pretty inefficient as it serves as a prototype.
def expand(arg_char, config):
    return naive_expand(arg_char, config)



## This function compiles an arg char by recursing if it contains quotes or command substitution.
##
## It is currently being extended to also expand any arguments that are safe to expand. 
def compile_arg_char(original_arg_char, fileIdGen, config):
    ## Check if the arg char can be expanded and if so do it.
    arg_char = copy.deepcopy(original_arg_char)
    if(safe_to_expand(arg_char)):
        arg_char = expand(arg_char, config)

    ## Compile the arg char
    key, val = get_kv(arg_char)
    if (key in ['C',   # Single character
                'E']): # Escape
        return arg_char
    elif (key == 'B'):
        ## TODO: I probably have to redirect the input of the compiled
        ##       node (IR) to be closed, and the output to be
        ##       redirected to some file that we will use to write to
        ##       the command argument to complete the command
        ##       substitution.
        compiled_node = compile_node(val, fileIdGen, config)
        return [key, compiled_node]
    elif (key == 'Q'):
        compiled_val = compile_command_argument(val, fileIdGen, config)
        return [key, compiled_val]
    else:
        ## TODO: Complete this
        return arg_char

def compile_command_argument(argument, fileIdGen, config):
    compiled_argument = [compile_arg_char(char, fileIdGen, config) for char in argument]
    return compiled_argument

def compile_command_arguments(arguments, fileIdGen, config):
    compiled_arguments = [compile_command_argument(arg, fileIdGen, config) for arg in arguments]
    return compiled_arguments

## Compiles the value assigned to a variable using the command argument rules.
## TODO: Is that the correct way to handle them?
def compile_assignments(assignments, fileIdGen, config):
    compiled_assignments = [[assignment[0], compile_command_argument(assignment[1], fileIdGen, config)]
                            for assignment in assignments]
    return compiled_assignments

def compile_redirection(redirection, fileIdGen, config):
    redir_type = redirection[0]
    redir_subtype = redirection[1][0]
    stream_id = redirection[1][1]
    file_arg = compile_command_argument(redirection[1][2], fileIdGen, config)
    return [redir_type, [redir_subtype, stream_id, file_arg]]

def compile_redirections(redirections, fileIdGen, config):
    compiled_redirections = [compile_redirection(redirection, fileIdGen, config)
                             for redirection in redirections]
    return compiled_redirections

##
## Preprocessing
##

## The preprocessing pass replaces all _candidate_ dataflow regions with
## calls to PaSh's runtime to let it establish if they are actually dataflow
## regions. The pass serializes all candidate dataflow regions:
## - A list of ASTs if at the top level or
## - an AST subtree if at a lower level
##
## The PaSh runtime then deserializes them, compiles them (if safe) and optimizes them.

## Replace candidate dataflow AST regions with calls to PaSh's runtime.
def replace_ast_regions(ast_objects, irFileGen, config):
    ## TODO: Copy the structure from compile_asts
    ##       since we want to merge multiple asts 
    ##       into one if they contain dataflow regions.

    ## TODO: Follow exactly the checks that are done in compile_asts
    ##       without actually compiling the asts to IRs.
    preprocessed_asts = []
    candidate_dataflow_region = []
    for i, ast_object in enumerate(ast_objects):
        # print("Preprocessing AST {}".format(i))
        # print(ast_object)

        ## NOTE: This could also replace all ASTs with calls to PaSh runtime.
        ##       There are a coupld issues with that:
        ##       1. We would have to make sure that state changes from PaSh runtime
        ##          affect the current shell. However, this probably has to be solved 
        ##          anyway except if we can *ensure* that no state changes can happen
        ##          in replaced parts.
        ##       2. Performance issues. Performance would be bad.

        ## Goals: This transformation can approximate in several directions.
        ##        1. Not replacing a candidate dataflow region.
        ##        2. Replacing a too large candidate region 
        ##           (making expansion not happen as late as possible)
        ##        3. Not replacing a maximal dataflow region, 
        ##           e.g. splitting a big one into two.
        ##        4. Replacing sections that are *certainly* not dataflow regions.
        ##           (This can only lead to performance issues.)
        ##
        ##        Which of the above can we hope to be precise with?
        ##        Can we have proofs indicating that we are not approximating those? 
        
        ## Preprocess ast by replacing subtrees with calls to runtime.
        ## - If the whole AST needs to be replaced (e.g. if it is a pipeline)
        ##   then the second output is true.
        ## - If the next AST needs to be replaced too (e.g. if the current one is a background)
        ##   then the third output is true
        output = preprocess_node(ast_object, irFileGen, config)
        preprocessed_ast, should_replace_whole_ast, is_non_maximal = output
        ## If the dataflow region is not maximal then it implies that the whole
        ## AST should be replaced.
        assert(not is_non_maximal or should_replace_whole_ast)

        ## If it isn't maximal then we just add it to the candidate
        if(is_non_maximal):
            candidate_dataflow_region.append(preprocessed_ast)
        else:
            ## If the current candidate dataflow region is non-empty
            ## it means that the previous AST was in the background so
            ## the current one has to be included in the process no matter what    
            if (len(candidate_dataflow_region) > 0):
                candidate_dataflow_region.append(preprocessed_ast)
                ## Since the current one is maximal (or not wholy replaced) 
                ## we close the candidate.
                replaced_ast = replace_df_region(candidate_dataflow_region, irFileGen, config)
                candidate_dataflow_region = []
                preprocessed_asts.append(replaced_ast)
            else:
                if(should_replace_whole_ast):
                    replaced_ast = replace_df_region([preprocessed_ast], irFileGen, config)
                    preprocessed_asts.append(replaced_ast)
                else:
                    preprocessed_asts.append(preprocessed_ast)

    ## Close the final dataflow region
    if(len(candidate_dataflow_region) > 0):
        replaced_ast = replace_df_region(candidate_dataflow_region, irFileGen, config)
        candidate_dataflow_region = []
        preprocessed_asts.append(replaced_ast)
    
    return preprocessed_asts

def preprocess_node(ast_object, irFileGen, config):
    global preprocess_cases
    return ast_match_untyped(ast_object, preprocess_cases, irFileGen, config)

## This preprocesses the AST node and also replaces it if it needs replacement .
## It is called by constructs that cannot be included in a dataflow region.
def preprocess_close_node(ast_object, irFileGen, config):
    output = preprocess_node(ast_object, irFileGen, config)
    preprocessed_ast, should_replace_whole_ast, _is_non_maximal = output
    if(should_replace_whole_ast):
        ## TODO: Maybe the first argument has to be a singular list?
        final_ast = replace_df_region(preprocessed_ast, irFileGen, config)
    else:
        final_ast = preprocessed_ast
    return final_ast

def preprocess_node_pipe(ast_node, _irFileGen, _config):
    ## A pipeline is *always* a candidate dataflow region.
    ## Q: Is that true?

    ## TODO: Preprocess the internals of the pipe to allow
    ##       for mutually recursive calls to PaSh.
    ##
    ##       For example, if a command in the pipe has a command substitution
    ##       in one of its arguments then we would like to call our runtime
    ##       there instead of 
    return ast_node, True, ast_node.is_background

## TODO: Complete this
def preprocess_node_command(ast_node, _irFileGen, _config):
    ## TODO: Preprocess the internals of the pipe to allow
    ##       for mutually recursive calls to PaSh.
    ##
    ##       For example, if a command in the pipe has a command substitution
    ##       in one of its arguments then we would like to call our runtime
    ##       there instead of 

    ## If there are no arguments, the command is just an
    ## assignment (Q: or just redirections?)
    if(len(ast_node.arguments) == 0):
        return ast_node, False, False
    
    ## This means we have a command. Commands are always candidate dataflow
    ## regions.
    return ast_node, True, False

## TODO: This is not efficient at all since it calls the PaSh runtime everytime the loop is entered.
##       We have to find a way to improve that.
def preprocess_node_for(ast_node, irFileGen, config):
    preprocessed_body = preprocess_close_node(ast_node.body, irFileGen, config)
    ## TODO: Could there be a problem with the in-place update
    ast_node.body = preprocessed_body
    return ast_node, False, False


## TODO: All the replace parts might need to be deleteD




## Replaces IR subtrees with a command that calls them (more
## precisely, a command that calls a python script to call them).
##
## Note: The traversal that replace_irs does, is exactly the same as
## the one that is done by compile_node. Both of these functions
## transform nodes of type t to something else.
##
## TODO: For now this just replaces the IRs starting from the ourside
## one first, but it should start from the bottom up to handle
## recursive IRs.
##
## TODO: Optimization: This pass could be merged with the previous
## translation one. I am not doing this now so that it remains
## cleaner. If we find out that we might benefit from this
## optimization, we can do it.
##
## TODO: Visual improvement: Can we abstract away the traverse of this
## ast in a mpa function, so that we don't rewrite all the cases?
def replace_irs(ast, irFileGen, config):

    if (isinstance(ast, IR)):
        replaced_ast = replace_df_region(ast, irFileGen, config)
    else:
        global ir_cases
        replaced_ast = ast_match(ast, ir_cases, irFileGen, config)

    return replaced_ast

## This function serializes a candidate df_region in a file, and in its place, 
## it adds a command that calls our distribution planner with the name of the
## saved file.
def replace_df_region(asts, irFileGen, config):
    ir_file_id = irFileGen.next_file_id()

    temp_ir_directory_prefix = config['distr_planner']['temp_ir_prefix']
    ## TODO: I probably have to generate random file names for the irs, so
    ## that multiple executions of PaSh don't interfere.
    ir_filename = ir_file_id.toFileName(temp_ir_directory_prefix)

    ## Serialize the node in a file
    with open(ir_filename, "wb") as ir_file:
        pickle.dump(asts, ir_file)

    ## Replace it with a command that calls the distribution
    ## planner with the name of the file.
    replaced_node = make_command(ir_filename)

    return replaced_node

## This function makes a command that calls the distribution planner
## together with the name of the file containing an IR. Then the
## distribution planner should read from this file and continue
## execution.
##
## (MAYBE) TODO: The way I did it, is by calling the parser once, and seeing
## what it returns. Maybe it would make sense to call the parser on
## the fly to have a cleaner implementation here?
def make_command(ir_filename):

    ## TODO: Do we need to do anything with the line_number? If so, make
    ## sure that I keep it in the IR, so that I can find it.
    arguments = [string_to_argument(config.PYTHON_VERSION),
                 string_to_argument(config.PLANNER_EXECUTABLE),
                 string_to_argument(ir_filename)]
    ## Pass a relevant argument to the planner
    arguments += config.pass_common_arguments(config.pash_args)

    line_number = 0
    node = make_kv('Command', [line_number, [], arguments, []])
    return node


def replace_irs_and_or_semi(ast_node, irFileGen, config):
    r_left = replace_irs(ast_node.left_operand, irFileGen, config)
    r_right = replace_irs(ast_node.right_operand, irFileGen, config)
    return make_kv(ast_node.construct.value, [r_left, r_right])

def replace_irs_command(ast_node, irFileGen, config):
    ## TODO: I probably have to also handle the redir_list (applies to
    ## compile_node_command too)
    c = ast_node.construct.value
    line_number = ast_node.line_number
    assignments = ast_node.assignments
    args = ast_node.arguments
    redir_list = ast_node.redir_list

    replaced_assignments = replace_irs_assignments(assignments, irFileGen, config)
    if(len(args) == 0):
        ## Note: The flow here is similar to compile
        replaced_ast = make_kv(c, [line_number, replaced_assignments, args, redir_list])
    else:
        command_name = args[0]
        replaced_options = replace_irs_command_arguments(args[1:], irFileGen, config)
        replaced_args = [command_name] + replaced_options

        replaced_ast = make_kv(c, [line_number, replaced_assignments, replaced_args, redir_list])

    return replaced_ast

def replace_irs_for(ast_node, irFileGen, config):
    ## Question: Is it a problem if we replace the same IR? Could
    ## there be changes between loops that are not reflected?
    replaced_ast = make_kv(ast_node.construct.value,
                           [ast_node.line_number,
                            replace_irs_command_argument(ast_node.argument, irFileGen, config),
                            replace_irs(ast_node.body, irFileGen, config),
                            ast_node.variable])
    return replaced_ast


def replace_irs_arg_char(arg_char, irFileGen, config):
    key, val = get_kv(arg_char)
    if (key == 'C'):
        return arg_char
    elif (key == 'B'):
        replaced_node = replace_irs(val, irFileGen, config)
        return make_kv(key, replaced_node)
    elif (key == 'Q'):
        replaced_val = replace_irs_command_argument(val, irFileGen, config)
        return make_kv(key, replaced_val)
    else:
        ## TODO: Complete this (as we have to do with the compile)
        return arg_char

def replace_irs_command_argument(argument, irFileGen, config):
    replaced_argument = [replace_irs_arg_char(char, irFileGen, config) for char in argument]
    return replaced_argument

def replace_irs_command_arguments(arguments, irFileGen, config):
    replaced_arguments = [replace_irs_command_argument(arg, irFileGen, config) for arg in arguments]
    return replaced_arguments

## This is similar to compile_assignments
##
## TODO: Is that the correct way to handle them? Question also applied
## to compile_assignments.
def replace_irs_assignments(assignments, irFileGen, config):
    replaced_assignments = [[assignment[0], replace_irs_command_argument(assignment[1], irFileGen, config)]
                            for assignment in assignments]
    return replaced_assignments

##
## Pattern matching for the AST
##

def check_if_ast_is_supported(construct, arguments, **kwargs):
    return

def ast_match_untyped(untyped_ast_object, cases, fileIdGen, config):
    ## TODO: This should construct the complete AstNode object (not just the surface level)
    ast_node = AstNode(untyped_ast_object)
    if ast_node.construct is AstNodeConstructor.PIPE:
        ast_node.check(children_count = lambda : len(ast_node.items) >= 2)
    return ast_match(ast_node, cases, fileIdGen, config)

def ast_match(ast_node, cases, fileIdGen, config):
    ## TODO: Remove that once `ast_match_untyped` is fixed to
    ##       construct the whole AstNode object.
    if(not isinstance(ast_node, AstNode)):
        return ast_match_untyped(ast_node, cases, fileIdGen, config)

    return cases[ast_node.construct.value](fileIdGen, config)(ast_node)
