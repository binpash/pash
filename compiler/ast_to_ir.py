from ir import *
from union_find import *
from definitions.ast_node import *
from definitions.ast_node_c import *
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
        if (command.category in ["stateless", "pure"]):
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


def compile_arg_char(arg_char, fileIdGen, config):
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
        replaced_ast = replace_ir(ast, irFileGen, config)
    else:
        global ir_cases
        replaced_ast = ast_match(ast, ir_cases, irFileGen, config)

    return replaced_ast

## This function serializes an IR in a file, and in its place, it adds
## a command that calls our distribution planner with the name of the
## saved file.
def replace_ir(ast_node, irFileGen, config):
    ir_file_id = irFileGen.next_file_id()

    temp_ir_directory_prefix = config['distr_planner']['temp_ir_prefix']
    ## TODO: I probably have to generate random file names for the irs, so
    ## that multiple executions of dish don't interfere.
    ir_filename = ir_file_id.toFileName(temp_ir_directory_prefix)

    ## Serialize the node in a file
    with open(ir_filename, "wb") as ir_file:
        pickle.dump(ast_node, ir_file)

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
    arguments += config.pass_common_arguments(config.dish_args)

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

def ast_match(ast_object, cases, fileIdGen, config):
    ast_node = AstNode(ast_object)
    if ast_node.construct is AstNodeConstructor.PIPE:
        ast_node.check(children_count = lambda : len(ast_node.items) >= 2)

    return cases[ast_node.construct.value](fileIdGen, config)(ast_node)
