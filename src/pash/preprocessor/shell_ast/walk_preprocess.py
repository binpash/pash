"""
Generic preprocessing walker for shell ASTs.

This module provides a pattern-matching based walker that separates
traversal mechanics from preprocessing policy.
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Callable, Any, TYPE_CHECKING

from shasta.ast_node import (
    AstNode,
    PipeNode,
    CommandNode,
    BackgroundNode,
    RedirNode,
    SubshellNode,
    SemiNode,
    AndNode,
    OrNode,
    NotNode,
    IfNode,
    ForNode,
    WhileNode,
    CaseNode,
    DefunNode,
    ArithNode,
    CondNode,
    SelectNode,
    ArithForNode,
    CoprocNode,
    TimeNode,
    GroupNode,
)

if TYPE_CHECKING:
    from shell_ast.transformation_options import AbstractTransformationState

from shell_ast.ast_util import PreprocessedAST


@dataclass
class PreprocessContext:
    """Context threaded through preprocessing traversal."""

    trans_options: Any  # AbstractTransformationState
    last_object: bool = False


@dataclass
class NodeResult:
    """Result from processing a single node."""

    ast: AstNode
    replace_whole: bool = False
    non_maximal: bool = False
    something_replaced: bool = False

    def to_preprocessed_ast(self, last_ast: bool) -> PreprocessedAST:
        """Convert to PreprocessedAST for API compatibility."""
        return PreprocessedAST(
            ast=self.ast,
            replace_whole=self.replace_whole,
            non_maximal=self.non_maximal,
            something_replaced=self.something_replaced,
            last_ast=last_ast,
        )


# Type alias for custom handlers
# Handler signature: (node, ctx, walker) -> NodeResult
NodeHandler = Callable[[AstNode, PreprocessContext, "WalkPreprocess"], NodeResult]


class WalkPreprocess:
    """
    Pattern-matching based preprocessing walker.

    Separates traversal mechanics from preprocessing policy through
    a handler-based architecture.
    """

    def __init__(self, handlers: dict[str, NodeHandler] | None = None):
        """
        Initialize the walker with optional custom handlers.

        Args:
            handlers: Dict mapping node type names (lowercase) to handler functions.
                     Handlers receive (node, ctx, walker) and return NodeResult.
        """
        self._handlers = handlers or {}

    def walk(self, node: AstNode, ctx: PreprocessContext) -> PreprocessedAST:
        """
        Walk and preprocess an AST node.

        Args:
            node: The node to preprocess
            ctx: The preprocessing context

        Returns:
            PreprocessedAST with transformed node and metadata
        """
        result = self._dispatch(node, ctx)
        return result.to_preprocessed_ast(ctx.last_object)

    def walk_close(
        self, node: AstNode, ctx: PreprocessContext
    ) -> tuple[AstNode, bool]:
        """
        Walk a node with "close" semantics - preprocess and optionally replace.

        This is used for children that cannot be part of the parent's
        dataflow region (e.g., children of control flow constructs).

        Args:
            node: The child node to preprocess
            ctx: The preprocessing context

        Returns:
            Tuple of (transformed_ast, something_replaced)
        """
        preprocessed = self.walk(node, ctx)

        if preprocessed.should_replace_whole_ast():
            final_ast = ctx.trans_options.replace_df_region(
                asts=[preprocessed.ast], disable_parallel_pipelines=ctx.last_object
            )
            return final_ast, True
        else:
            return preprocessed.ast, preprocessed.will_anything_be_replaced()

    def _dispatch(self, node: AstNode, ctx: PreprocessContext) -> NodeResult:
        """Dispatch to appropriate handler based on node type."""
        node_name = type(node).NodeName.lower()

        # Check for custom handler first
        if node_name in self._handlers:
            return self._handlers[node_name](node, ctx, self)

        # Fall back to default pattern-matching walker
        return self._default_walk(node, ctx)

    def _default_walk(self, node: AstNode, ctx: PreprocessContext) -> NodeResult:
        """Default walking behavior using pattern matching."""
        match node:
            # === Leaf replacement nodes ===
            # Note: When replace_whole=True, something_replaced must also be True
            # (this is asserted in ast_to_ast.py)
            case PipeNode():
                return NodeResult(
                    ast=node,
                    replace_whole=True,
                    non_maximal=node.is_background,
                    something_replaced=True,
                )

            case CommandNode():
                if len(node.arguments) == 0:
                    # Just an assignment, not a candidate
                    return NodeResult(ast=node, something_replaced=False)
                return NodeResult(
                    ast=node, replace_whole=True, something_replaced=True
                )

            case BackgroundNode():
                return NodeResult(
                    ast=node,
                    replace_whole=True,
                    non_maximal=True,
                    something_replaced=True,
                )

            # === Binary operators with close-node semantics ===
            case SemiNode():
                return self._walk_binary_close(
                    node, ctx, "left_operand", "right_operand"
                )

            case AndNode():
                return self._walk_binary_close(
                    node, ctx, "left_operand", "right_operand"
                )

            case OrNode():
                return self._walk_binary_close(
                    node, ctx, "left_operand", "right_operand"
                )

            # === Single-child close nodes ===
            case RedirNode():
                return self._walk_single_close(node, ctx, "node")

            case SubshellNode():
                return self._walk_single_close(node, ctx, "body")

            case NotNode():
                return self._walk_single_close(node, ctx, "body")

            case GroupNode():
                return self._walk_single_close(node, ctx, "body")

            case CoprocNode():
                return self._walk_single_close(node, ctx, "body")

            case TimeNode():
                return self._walk_single_close(node, ctx, "command")

            case SelectNode():
                return self._walk_single_close(node, ctx, "body")

            case ArithForNode():
                return self._walk_single_close(node, ctx, "action")

            # === Control flow with multiple children ===
            case WhileNode():
                return self._walk_while(node, ctx)

            case ForNode():
                return self._walk_for(node, ctx)

            case IfNode():
                return self._walk_if(node, ctx)

            case CaseNode():
                return self._walk_case(node, ctx)

            case CondNode():
                return self._walk_cond(node, ctx)

            # === No-op nodes ===
            case DefunNode():
                return NodeResult(ast=node, something_replaced=False)

            case ArithNode():
                return NodeResult(ast=node, something_replaced=False)

            case _:
                raise ValueError(f"Unknown node type: {type(node).NodeName}")

    # === Helper methods for common patterns ===

    def _walk_single_close(
        self, node: AstNode, ctx: PreprocessContext, child_attr: str
    ) -> NodeResult:
        """Walk a node with a single child using close-node semantics."""
        child = getattr(node, child_attr)
        new_child, replaced = self.walk_close(child, ctx)
        setattr(node, child_attr, new_child)
        return NodeResult(ast=node, something_replaced=replaced)

    def _walk_binary_close(
        self,
        node: AstNode,
        ctx: PreprocessContext,
        left_attr: str,
        right_attr: str,
    ) -> NodeResult:
        """Walk a binary node with close-node semantics on both children."""
        left = getattr(node, left_attr)
        right = getattr(node, right_attr)

        new_left, replaced_left = self.walk_close(left, ctx)
        new_right, replaced_right = self.walk_close(right, ctx)

        setattr(node, left_attr, new_left)
        setattr(node, right_attr, new_right)

        return NodeResult(ast=node, something_replaced=replaced_left or replaced_right)

    def _walk_while(self, node: WhileNode, ctx: PreprocessContext) -> NodeResult:
        """Walk a while node."""
        # Enter loop context
        ctx.trans_options.enter_loop()

        new_test, test_replaced = self.walk_close(node.test, ctx)
        new_body, body_replaced = self.walk_close(node.body, ctx)

        node.test = new_test
        node.body = new_body

        # Exit loop context
        ctx.trans_options.exit_loop()

        return NodeResult(
            ast=node, something_replaced=test_replaced or body_replaced
        )

    def _walk_for(self, node: ForNode, ctx: PreprocessContext) -> NodeResult:
        """Walk a for node (base behavior without loop tracking injection)."""
        # Enter loop context
        ctx.trans_options.enter_loop()

        new_body, replaced = self.walk_close(node.body, ctx)
        node.body = new_body

        # Exit loop context
        ctx.trans_options.exit_loop()

        return NodeResult(ast=node, something_replaced=replaced)

    def _walk_if(self, node: IfNode, ctx: PreprocessContext) -> NodeResult:
        """Walk an if node."""
        new_cond, cond_replaced = self.walk_close(node.cond, ctx)

        ctx.trans_options.enter_if()
        new_then, then_replaced = self.walk_close(node.then_b, ctx)

        if node.else_b is not None:
            ctx.trans_options.enter_else()
            new_else, else_replaced = self.walk_close(node.else_b, ctx)
        else:
            new_else, else_replaced = None, False

        ctx.trans_options.exit_if()

        node.cond = new_cond
        node.then_b = new_then
        node.else_b = new_else

        return NodeResult(
            ast=node,
            something_replaced=cond_replaced or then_replaced or else_replaced,
        )

    def _walk_case(self, node: CaseNode, ctx: PreprocessContext) -> NodeResult:
        """Walk a case node."""
        any_replaced = False
        new_cases = []

        for case in node.cases:
            if case.get("cbody") is not None:
                new_body, replaced = self.walk_close(case["cbody"], ctx)
                case["cbody"] = new_body
                any_replaced = any_replaced or replaced
            new_cases.append(case)

        node.cases = new_cases
        return NodeResult(ast=node, something_replaced=any_replaced)

    def _walk_cond(self, node: CondNode, ctx: PreprocessContext) -> NodeResult:
        """Walk a cond node (bash [[ ]] expression)."""
        replaced_left = False
        replaced_right = False

        if node.left is not None:
            node.left, replaced_left = self.walk_close(node.left, ctx)
        if node.right is not None:
            node.right, replaced_right = self.walk_close(node.right, ctx)

        return NodeResult(
            ast=node, something_replaced=replaced_left or replaced_right
        )
