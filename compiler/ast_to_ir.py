import subprocess

from shasta.ast_node import *
from sh_expand.expand import expand_command, ExpansionState

from shell_ast.ast_util import *
from ir import *
from util import *
from parse import from_ast_objects_to_shell

from custom_error import *

## TODO: Separate the ir stuff to the bare minimum and
##       try to move this to the shell_ast folder.

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
    "Pipe": (
        lambda fileIdGen, config: lambda ast_node: compile_node_pipe(
            ast_node, fileIdGen, config
        )
    ),
    "Command": (
        lambda fileIdGen, config: lambda ast_node: compile_node_command(
            ast_node, fileIdGen, config
        )
    ),
    "And": (
        lambda fileIdGen, config: lambda ast_node: compile_node_and_or_semi(
            ast_node, fileIdGen, config
        )
    ),
    "Or": (
        lambda fileIdGen, config: lambda ast_node: compile_node_and_or_semi(
            ast_node, fileIdGen, config
        )
    ),
    "Semi": (
        lambda fileIdGen, config: lambda ast_node: compile_node_and_or_semi(
            ast_node, fileIdGen, config
        )
    ),
    "Redir": (
        lambda fileIdGen, config: lambda ast_node: compile_node_redir_subshell(
            ast_node, fileIdGen, config
        )
    ),
    "Subshell": (
        lambda fileIdGen, config: lambda ast_node: compile_node_redir_subshell(
            ast_node, fileIdGen, config
        )
    ),
    "Background": (
        lambda fileIdGen, config: lambda ast_node: compile_node_background(
            ast_node, fileIdGen, config
        )
    ),
    "For": (
        lambda fileIdGen, config: lambda ast_node: compile_node_for(
            ast_node, fileIdGen, config
        )
    ),
}


def compile_asts(ast_objects: "list[AstNode]", fileIdGen, config):
    compiled_asts = []
    acc_ir = None
    for i, ast_object in enumerate(ast_objects):
        # log("Compiling AST {}".format(i))
        # log(ast_object)
        assert isinstance(ast_object, AstNode)

        ## Compile subtrees of the AST to out intermediate representation
        ## KK 2023-05-25: Would we ever want to pass this state to the expansion
        ##                of the next object? I don't think so.
        exp_state = ExpansionState(config["shell_variables"])
        expanded_ast = expand_command(ast_object, exp_state)
        # log("Expanded:", expanded_ast)
        compiled_ast = compile_node(expanded_ast, fileIdGen, config)

        # log("Compiled AST:")
        # log(compiled_ast)

        ## If the accumulator contains an IR (meaning that the
        ## previous commands where run in background), union it with
        ## the current returned ast.
        if not acc_ir is None:
            if isinstance(compiled_ast, IR):
                acc_ir.background_union(compiled_ast)
            else:
                ## TODO: Make this union the compiled_ast with the
                ## accumulated IR, since the user wanted to run these
                ## commands in parallel (Is that correct?)
                # acc_ir.background_union(IR([compiled_ast]))
                compiled_asts.append(acc_ir)
                acc_it = None
                compiled_asts.append(compiled_ast)

            ## If the current compiled ast not in background (and so
            ## the union isn't in background too), stop accumulating
            if not acc_ir is None and not acc_ir.is_in_background():
                compiled_asts.append(acc_ir)
                acc_ir = None
        else:
            ## If the compiled ast is in background, start
            ## accumulating it
            if isinstance(compiled_ast, IR) and compiled_ast.is_in_background():
                acc_ir = compiled_ast
            else:
                compiled_asts.append(compiled_ast)

    ## The final accumulator
    if not acc_ir is None:
        compiled_asts.append(acc_ir)

    return compiled_asts


def compile_node(ast_object, fileIdGen, config):
    global compile_cases
    return ast_match(ast_object, compile_cases, fileIdGen, config)


def compile_node_pipe(ast_node, fileIdGen, config):
    compiled_pipe_nodes = combine_pipe(
        [compile_node(pipe_item, fileIdGen, config) for pipe_item in ast_node.items]
    )

    ## Note: When calling combine_pipe_nodes (which
    ##       optimistically distributes all the children of a
    ##       pipeline) the compiled_pipe_nodes should always
    ##       be one IR
    compiled_ir = compiled_pipe_nodes[0]
    ## Save the old ast for the end-to-end prototype
    old_untyped_ast_node = ast_node.json()
    compiled_ir.set_ast(old_untyped_ast_node)
    ## Set the IR background so that it can be parallelized with
    ## the next command if the pipeline was run in background
    compiled_ir.set_background(ast_node.is_background)
    compiled_ast = compiled_ir
    return compiled_ast


## This combines all the children of the Pipeline to an IR.
def combine_pipe(ast_nodes):
    ## Initialize the IR with the first node in the Pipe
    if isinstance(ast_nodes[0], IR):
        combined_nodes = ast_nodes[0]
    else:
        ## If any part of the pipe is not an IR, the compilation must fail.
        log("Node: {} is not pure".format(ast_nodes[0]))
        raise UnparallelizableError("Node: {} is not a pure node in pipe".format(ast_nodes[0]))

    ## Combine the rest of the nodes
    for ast_node in ast_nodes[1:]:
        if isinstance(ast_node, IR):
            combined_nodes.pipe_append(ast_node)
        else:
            ## If any part of the pipe is not an IR, the compilation must fail.
            log("Node: {} is not pure".format(ast_nodes))
            raise UnparallelizableError("This specific node: {} is not a pure node in pipe".format(ast_node))

    return [combined_nodes]


def compile_node_command(ast_node, fileIdGen, config):
    ## Compile assignments and redirection list
    compiled_assignments = compile_assignments(ast_node.assignments, fileIdGen, config)
    compiled_redirections = compile_redirections(ast_node.redir_list, fileIdGen, config)

    ## This should never be possible since the preprocessor
    ##  wouldn't replace a call without arguments (simple assignment).
    assert len(ast_node.arguments) > 0

    arguments = ast_node.arguments
    command_name = arguments[0]
    options = compile_command_arguments(arguments[1:], fileIdGen, config)

    try:
        ## If the command is not compileable to a DFG the following call will fail
        ir = compile_command_to_DFG(
            fileIdGen, command_name, options, redirections=compiled_redirections
        )
        compiled_ast = ir
    except ValueError as err:
        log("Command not compiled to DFG:", err)
        ## TODO: Maybe we want to fail here instead of waiting for later?
        ##       Is there any case where a non-compiled command is fine?
        # log(traceback.format_exc())
        compiled_arguments = compile_command_arguments(arguments, fileIdGen, config)
        compiled_ast = make_kv(
            type(ast_node).NodeName,
            [
                ast_node.line_number,
                compiled_assignments,
                compiled_arguments,
                compiled_redirections,
            ],
        )

    return compiled_ast


def compile_node_and_or_semi(ast_node, fileIdGen, config):
    compiled_ast = make_kv(
        type(ast_node).NodeName,
        [
            compile_node(ast_node.left_operand, fileIdGen, config),
            compile_node(ast_node.right_operand, fileIdGen, config),
        ],
    )
    return compiled_ast


def compile_node_redir_subshell(ast_node, fileIdGen, config):
    compiled_node = compile_node(ast_node.node, fileIdGen, config)

    if isinstance(compiled_node, IR):
        ## TODO: I should use the redir list to redirect the files of
        ##       the IR accordingly
        compiled_ast = compiled_node
    else:
        compiled_ast = make_kv(
            type(ast_node).NodeName,
            [ast_node.line_number, compiled_node, ast_node.redir_list],
        )

    return compiled_ast


def compile_node_background(ast_node, fileIdGen, config):
    compiled_node = compile_node(ast_node.node, fileIdGen, config)

    ## TODO: I should use the redir list to redirect the files of
    ##       the IR accordingly
    if isinstance(compiled_node, IR):
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
        ## TODO: We should not compile the ast here
        compiled_ast = ast_node
        compiled_ast.node = compiled_node

    return compiled_ast


def compile_node_for(ast_node, fileIdGen, config):
    ## TODO: Investigate what kind of check could we do to make a for
    ## loop parallel
    compiled_ast = make_kv(
        type(ast_node).NodeName,
        [
            ast_node.line_number,
            compile_command_argument(ast_node.argument, fileIdGen, config),
            compile_node(ast_node.body, fileIdGen, config),
            ast_node.variable,
        ],
    )
    return compiled_ast


## This function checks if we should expand an arg_char
##
## It has a dual purpose:
## 1. First, it checks whether we need
##    to pay the overhead of expanding (by checking that there is something to expand)
##    like a variable
## 2. Second it raises an error if we cannot expand an argument.
def should_expand_arg_char(arg_char):
    key, val = get_kv(arg_char)
    if key in ["V"]:  # Variable
        return True
    elif key == "Q":
        return should_expand_argument(val)
    elif key == "B":
        log("Cannot expand:", arg_char)
        raise NotImplementedError()
    return False


def should_expand_argument(argument):
    return any([should_expand_arg_char(arg_char) for arg_char in argument])


## TODO: Move this function somewhere more general
def execute_shell_asts(asts):
    output_script = from_ast_objects_to_shell(asts)
    # log(output_script)
    exec_obj = subprocess.run(
        ["/usr/bin/env", "bash"],
        input=output_script,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
    )
    exec_obj.check_returncode()
    # log(exec_obj.stdout)
    return exec_obj.stdout


## TODO: Properly parse the output of the shell script
def parse_string_to_arguments(arg_char_string):
    # log(arg_char_string)
    return string_to_arguments(arg_char_string)


## TODO: Use "pash_input_args" when expanding in place of normal arguments.
def naive_expand(argument, config):
    ## config contains a dictionary with:
    ##  - all variables, their types, and values in 'shell_variables'
    ##  - the name of a file that contains them in 'shell_variables_file_path'
    # log(config['shell_variables'])
    # log(config['shell_variables_file_path'])

    ## Create an AST node that "echo"s the argument
    echo_asts = make_echo_ast(argument, config["shell_variables_file_path"])

    ## Execute the echo AST by unparsing it to shell
    ## and calling bash
    expanded_string = execute_shell_asts(echo_asts)

    log("Argument:", argument, "was expanded to:", expanded_string)

    ## Parse the expanded string back to an arg_char
    expanded_arguments = parse_string_to_arguments(expanded_string)

    ## TODO: Handle any errors
    # log(expanded_arg_char)
    return expanded_arguments


## This function expands an arg_char.
## At the moment it is pretty inefficient as it serves as a prototype.
##
## TODO: At the moment this has the issue that a command that has the words which we want to expand
##       might have assignments of its own, therefore requiring that we use them to properly expand.
def expand_command_argument(argument, config):
    new_arguments = [argument]
    if should_expand_argument(argument):
        new_arguments = naive_expand(argument, config)
    return new_arguments


## This function compiles an arg char by recursing if it contains quotes or command substitution.
##
## It is currently being extended to also expand any arguments that are safe to expand.
def compile_arg_char(arg_char: ArgChar, fileIdGen, config):
    ## Compile the arg char
    if isinstance(arg_char, CArgChar) or isinstance(arg_char, EArgChar):
        # Single character or escape
        return arg_char
    elif isinstance(arg_char, BArgChar):
        ## TODO: I probably have to redirect the input of the compiled
        ##       node (IR) to be closed, and the output to be
        ##       redirected to some file that we will use to write to
        ##       the command argument to complete the command
        ##       substitution.
        arg_char.node = compile_node(arg_char.node, fileIdGen, config)
        return arg_char
    elif isinstance(arg_char, QArgChar):
        arg_char.arg = compile_command_argument(arg_char.arg, fileIdGen, config)
        return arg_char
    else:
        log(f"Unknown arg_char: {arg_char}")
        ## TODO: Complete this
        return arg_char


def compile_command_argument(argument, fileIdGen, config):
    compiled_argument = [compile_arg_char(char, fileIdGen, config) for char in argument]
    return compiled_argument


def compile_command_arguments(arguments, fileIdGen, config):
    compiled_arguments = [
        compile_command_argument(arg, fileIdGen, config) for arg in arguments
    ]
    return compiled_arguments


## Compiles the value assigned to a variable using the command argument rules.
## TODO: Is that the correct way to handle them?
def compile_assignments(assignments, fileIdGen, config):
    compiled_assignments = [
        [assignment[0], compile_command_argument(assignment[1], fileIdGen, config)]
        for assignment in assignments
    ]
    return compiled_assignments


def compile_redirection(redirection, fileIdGen, config):
    file_arg = compile_command_argument(redirection.arg, fileIdGen, config)
    redirection.arg = file_arg
    return redirection


def compile_redirections(redirections, fileIdGen, config):
    compiled_redirections = [
        compile_redirection(redirection, fileIdGen, config)
        for redirection in redirections
    ]
    return compiled_redirections
