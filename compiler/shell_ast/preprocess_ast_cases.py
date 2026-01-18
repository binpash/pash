"""
AST preprocessing for PaSh.

This module uses the generic walk_preprocess infrastructure with
custom handlers for PaSh-specific behavior.
"""

from shell_ast.walk_preprocess import WalkPreprocess, PreprocessContext
from shell_ast.handlers.loop_tracking import for_node_with_loop_tracking
from shell_ast.ast_util import PreprocessedAST
from shell_ast.transformation_options import AbstractTransformationState
from shasta.ast_node import AstNode


def create_pash_walker() -> WalkPreprocess:
    """
    Create a preprocessing walker configured for PaSh.

    Returns:
        WalkPreprocess instance with PaSh-specific handlers
    """
    handlers = {
        "for": for_node_with_loop_tracking,
    }
    return WalkPreprocess(handlers=handlers)


# Module-level walker instance
_pash_walker = create_pash_walker()


def preprocess_node(
    ast_node: AstNode,
    trans_options: AbstractTransformationState,
    last_object: bool,
) -> PreprocessedAST:
    """
    Preprocesses an AstNode. Given an AstNode of any type, it will appropriately
    dispatch a preprocessor for the specific node type.

    This is the main entry point for preprocessing, maintaining backward
    compatibility with the existing API.

    Parameters:
        ast_node (AstNode): The AstNode to parse
        trans_options (AbstractTransformationState):
            A concrete transformation state instance corresponding to the output target
        last_object (bool): Flag for whether this is the last AstNode

    Returns:
        PreprocessedAst: the preprocessed version of the original AstNode
    """
    ctx = PreprocessContext(trans_options=trans_options, last_object=last_object)
    return _pash_walker.walk(ast_node, ctx)


def preprocess_close_node(
    ast_node: AstNode,
    trans_options: AbstractTransformationState,
    last_object: bool = False,
):
    """
    Preprocess an AST node and replace if needed.

    This preprocesses the AST node and also replaces it if it needs replacement.
    It is called by constructs that cannot be included in a dataflow region.

    Parameters:
        ast_node (AstNode): The AstNode to preprocess
        trans_options (AbstractTransformationState):
            A concrete transformation state instance
        last_object (bool): Flag for whether this is the last AstNode

    Returns:
        Tuple of (final_ast, something_replaced)
    """
    ctx = PreprocessContext(trans_options=trans_options, last_object=last_object)
    return _pash_walker.walk_close(ast_node, ctx)
