import os
import subprocess

from shasta.ast_node import (
    AstNode,
    ArgChar,
    CArgChar,
    EArgChar,
    BArgChar,
    QArgChar,
    PipeNode,
    CommandNode,
    BackgroundNode,
)
from shasta.ast_walker import walk_ast_node
from sh_expand.expand import expand_command, ExpansionState

import config
from sh_expand.bash_expand import (
    BashExpansionState,
    expand_command as expand_command_bash,
)
from ir import IR, compile_command_to_DFG
from parse import from_ast_objects_to_shell
from util import (
    log,
    get_kv,
    make_kv,
    string_to_arguments,
    UnparallelizableError,
)


BASH_MODE = os.environ.get("pash_shell") == "bash"
if BASH_MODE:
    # doesn't need to be kept alive accross invocations
    # but should be faster to avoid creating new processes

    # too verbose to keep on with -d 1 (TODO: Get -d 2 working here)
    BASH_EXP_STATE = BashExpansionState(temp_dir=config.PASH_TMP_PREFIX, debug=False)

##
## Compile AST -> Extended AST with IRs
##
## The preprocessor sends only dataflow region nodes to the compiler:
## - PipeNode
## - CommandNode (with arguments)
## - BackgroundNode
##
## All other node types (control flow, wrappers, etc.) stay in the
## preprocessed script and are executed by bash directly.
##


def compile_asts(ast_objects: "list[AstNode]", fileIdGen, config):
    compiled_asts = []
    acc_ir = None
    for i, ast_object in enumerate(ast_objects):
        assert isinstance(ast_object, AstNode)

        ## Compile subtrees of the AST to our intermediate representation
        if BASH_MODE:
            if not BASH_EXP_STATE.is_open:
                BASH_EXP_STATE.open()
            expanded_ast = expand_command_bash(
                ast_object,
                BASH_EXP_STATE,
                config["shell_variables_file_path"],
                config["shell_variables"],
            )
        else:
            exp_state = ExpansionState(config["shell_variables"])
            expanded_ast = expand_command(ast_object, exp_state)
        compiled_ast = compile_node(expanded_ast, fileIdGen, config)

        ## If the accumulator contains an IR (meaning that the
        ## previous commands where run in background), union it with
        ## the current returned ast.
        if not acc_ir is None:
            if isinstance(compiled_ast, IR):
                acc_ir.background_union(compiled_ast)
            else:
                ## shouldn't happen since compile_node should have already
                ## raised this error
                raise UnparallelizableError(f"Node: {compiled_ast} is not pure")

            ## If the current compiled ast not in background (and so
            ## the union isn't in background too), stop accumulating
            if not acc_ir.is_in_background():
                compiled_asts.append(acc_ir)
                acc_ir = None
        else:
            ## If the compiled ast is in background, start
            ## accumulating it
            if compiled_ast.is_in_background():
                acc_ir = compiled_ast
            else:
                compiled_asts.append(compiled_ast)

    ## The final accumulator
    if not acc_ir is None:
        compiled_asts.append(acc_ir)

    return compiled_asts


def make_compile_transform(fileIdGen, config):
    """
    Create a compile transformation function for use with walk_ast_node.

    Only handles node types that actually reach the compiler from the preprocessor:
    - PipeNode, CommandNode, BackgroundNode (dataflow region nodes)
    - BArgChar, QArgChar (for command substitution and quoted strings)
    """

    def compile_transform(node):
        match node:
            case PipeNode():
                # Compile pipe children and combine into single IR
                compiled_items = [compile_node(item, fileIdGen, config) for item in node.items]
                compiled_ir = combine_pipe(compiled_items)
                compiled_ir.set_ast(node.json())
                compiled_ir.set_background(node.is_background)
                return compiled_ir

            case CommandNode():
                return compile_command_node(node, fileIdGen, config)

            case BackgroundNode():
                compiled_inner = compile_node(node.node, fileIdGen, config)
                if isinstance(compiled_inner, IR):
                    compiled_inner.set_background(True)
                    return compiled_inner
                else:
                    # Note: background nodes can be added in the distributed graph
                    # similarly to the children of pipelines.
                    node.node = compiled_inner
                    return node

            case BArgChar():
                # Command substitution - compile the inner node
                return BArgChar(node=compile_node(node.node, fileIdGen, config))

            case QArgChar():
                # Quoted string - compile the inner argument
                return QArgChar(arg=compile_command_argument(node.arg, fileIdGen, config))

            case _:
                # Let walker handle default recursion for other nodes
                return None

    return compile_transform


def compile_node(ast_object, fileIdGen, config) -> IR:
    """
    Compile an AST node to an IR or transformed AST.

    Uses walk_ast_node with a compile transformation function that handles
    the IR-specific logic for each node type.
    """
    transform = make_compile_transform(fileIdGen, config)
    return walk_ast_node(ast_object, replace=transform)


def combine_pipe(compiled_items):
    """
    Combine all compiled children of a Pipeline into a single IR.

    All children must be IRs; raises UnparallelizableError otherwise.
    """
    if not isinstance(compiled_items[0], IR):
        log("Node: {} is not pure".format(compiled_items[0]))
        raise UnparallelizableError("Node: {} is not a pure node in pipe".format(compiled_items[0]))

    combined = compiled_items[0]
    for item in compiled_items[1:]:
        if isinstance(item, IR):
            combined.pipe_append(item)
        else:
            log("Node: {} is not pure".format(item))
            raise UnparallelizableError("This specific node: {} is not a pure node in pipe".format(item))

    return combined


def compile_command_node(ast_node, fileIdGen, config):
    """
    Compile a CommandNode to IR or transformed AST.

    Tries to compile to a DFG; if that fails, returns the command
    wrapped in make_kv format.
    """
    # Compile assignments and redirection list
    compiled_assignments = compile_assignments(ast_node.assignments, fileIdGen, config)
    compiled_redirections = compile_redirections(ast_node.redir_list, fileIdGen, config)

    # This should never be possible since the preprocessor
    # wouldn't replace a call without arguments (simple assignment).
    assert len(ast_node.arguments) > 0

    arguments = ast_node.arguments
    command_name = arguments[0]
    options = compile_command_arguments(arguments[1:], fileIdGen, config)

    try:
        # If the command is not compileable to a DFG the following call will fail
        ir = compile_command_to_DFG(
            fileIdGen, command_name, options, redirections=compiled_redirections
        )
        return ir
    except ValueError as err:
        log("Command not compiled to DFG:", err)
        compiled_arguments = compile_command_arguments(arguments, fileIdGen, config)
        return make_kv(
            type(ast_node).NodeName,
            [
                ast_node.line_number,
                compiled_assignments,
                compiled_arguments,
                compiled_redirections,
            ],
        )


## This function checks if we should expand an arg_char
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


def execute_shell_asts(asts):
    output_script = from_ast_objects_to_shell(asts)
    exec_obj = subprocess.run(
        ["/usr/bin/env", "bash"],
        input=output_script,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
    )
    exec_obj.check_returncode()
    return exec_obj.stdout


def parse_string_to_arguments(arg_char_string):
    return string_to_arguments(arg_char_string)


def compile_arg_char(arg_char: ArgChar, fileIdGen, config):
    """Compile an arg char by recursing if it contains quotes or command substitution."""
    if isinstance(arg_char, CArgChar) or isinstance(arg_char, EArgChar):
        return arg_char
    elif isinstance(arg_char, BArgChar):
        arg_char.node = compile_node(arg_char.node, fileIdGen, config)
        return arg_char
    elif isinstance(arg_char, QArgChar):
        arg_char.arg = compile_command_argument(arg_char.arg, fileIdGen, config)
        return arg_char
    else:
        log(f"Unknown arg_char: {arg_char}")
        return arg_char


def compile_command_argument(argument, fileIdGen, config):
    return [compile_arg_char(char, fileIdGen, config) for char in argument]


def compile_command_arguments(arguments, fileIdGen, config):
    return [compile_command_argument(arg, fileIdGen, config) for arg in arguments]


def compile_assignment(assignment, fileIdGen, config):
    assignment.val = compile_command_argument(assignment.val, fileIdGen, config)
    return assignment


def compile_assignments(assignments, fileIdGen, config):
    return [compile_assignment(assignment, fileIdGen, config) for assignment in assignments]


def compile_redirection(redirection, fileIdGen, config):
    if redirection.arg is not None:
        redirection.arg = compile_command_argument(redirection.arg, fileIdGen, config)
    return redirection


def compile_redirections(redirections, fileIdGen, config):
    return [compile_redirection(redirection, fileIdGen, config) for redirection in redirections]
