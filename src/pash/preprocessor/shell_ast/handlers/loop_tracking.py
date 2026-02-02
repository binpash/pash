"""
Handler for ForNode that injects loop tracking code.

This is specific to the pash/compiler version which needs to
track loop iterations for speculative execution.
"""

from __future__ import annotations
import copy
from typing import TYPE_CHECKING

from shasta.ast_node import ForNode, make_typed_semi_sequence
from shasta.json_to_ast import to_ast_node
from shell_ast.ast_util import (
    make_export_var_constant_string,
    make_increment_var,
    export_pash_loop_iters_for_current_context,
)
from env_var_names import loop_iter_var

if TYPE_CHECKING:
    from shell_ast.walk_preprocess import PreprocessContext, NodeResult, WalkPreprocess


def for_node_with_loop_tracking(
    node: ForNode, ctx: "PreprocessContext", walker: "WalkPreprocess"
) -> "NodeResult":
    """
    Handler for ForNode that injects loop tracking code.

    This handler:
    1. Enters a new loop context
    2. Preprocesses the loop body
    3. Injects loop iteration counter initialization
    4. Injects loop iteration counter increment in the body
    5. Exports loop iteration context for the runtime
    6. Exits the loop context

    Args:
        node: The ForNode to process
        ctx: The preprocessing context
        walker: The walker instance for recursive calls

    Returns:
        NodeResult with the transformed node and metadata
    """
    # Import here to avoid circular dependency
    from shell_ast.walk_preprocess import NodeResult

    # Enter loop context
    loop_id = ctx.trans_options.enter_loop()

    # Preprocess the body using close-node semantics
    preprocessed_body, something_replaced = walker.walk_close(node.body, ctx)

    # Create loop tracking nodes
    var_name = loop_iter_var(loop_id)
    export_node = make_export_var_constant_string(var_name, "0")
    increment_node = make_increment_var(var_name)

    # Get all loop IDs for context export
    all_loop_ids = ctx.trans_options.get_current_loop_context()
    save_loop_iters_node = export_pash_loop_iters_for_current_context(all_loop_ids)

    # Modify the loop body to include tracking
    node.body = make_typed_semi_sequence(
        [
            to_ast_node(increment_node),
            to_ast_node(save_loop_iters_node),
            ## KK 2026-01-20 This deepcopy used to be there but led to serious
            ##               performance issues with some tests. 
            ##               I think it is not necessary but leaving it here in case we want to revert back.
            # copy.deepcopy(preprocessed_body),
            preprocessed_body,
        ]
    )

    # Exit loop context
    ctx.trans_options.exit_loop()

    # Reset loop iters after exiting
    out_of_loop_ids = ctx.trans_options.get_current_loop_context()
    reset_loop_iters_node = export_pash_loop_iters_for_current_context(out_of_loop_ids)

    # Wrap the entire for loop with initialization and reset
    new_node = make_typed_semi_sequence(
        [to_ast_node(export_node), node, to_ast_node(reset_loop_iters_node)]
    )

    return NodeResult(ast=new_node, something_replaced=something_replaced)
