"""
Preprocessing walker for shell ASTs.

This module provides a preprocessing visitor that builds on the generic
:class:`shasta.ast_walker.CommandVisitor` to add PaSh-specific
dataflow-region detection and replacement logic.
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import Any

from shasta.ast_node import (
    AstNode,
    PipeNode,
    CommandNode,
    BackgroundNode,
    CaseNode,
    DefunNode,
    ArithNode,
)
from shasta.ast_walker import CommandVisitor, command_child_attrs

from ast_util import PreprocessedAST


@dataclass
class PreprocessContext:
    """Context threaded through preprocessing traversal."""

    trans_options: Any  # TransformationState
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


class WalkPreprocess(CommandVisitor):
    """
    Preprocessing visitor for shell ASTs.

    Extends :class:`CommandVisitor` with preprocessing-specific behaviour:
    dataflow-region detection and close-node semantics.

    Most node types use the default :meth:`generic_visit` which walks
    all command children with close-node semantics.  Only nodes with
    truly specific behaviour (leaves, no-ops) override.
    """

    def __init__(self):
        self.ctx: PreprocessContext | None = None

    # === Public API ===

    def walk(self, node: AstNode, ctx: PreprocessContext) -> PreprocessedAST:
        """Walk and preprocess an AST node."""
        self.ctx = ctx
        result = self.visit(node)
        return result.to_preprocessed_ast(ctx.last_object)

    def walk_close(
        self, node: AstNode, ctx: PreprocessContext
    ) -> tuple[AstNode, bool]:
        """
        Walk a node with "close" semantics - preprocess and optionally replace.

        Used for children that cannot be part of the parent's dataflow region
        (e.g., children of control flow constructs).
        """
        preprocessed = self.walk(node, ctx)

        if preprocessed.should_replace_whole_ast():
            final_ast = ctx.trans_options.replace_df_region(
                asts=[preprocessed.ast], disable_parallel_pipelines=ctx.last_object
            )
            return final_ast, True
        else:
            return preprocessed.ast, preprocessed.will_anything_be_replaced()

    # === Default: walk all children with close semantics ===

    def generic_visit(self, node: AstNode) -> NodeResult:
        """Walk all command children with close-node semantics."""
        ctx = self.ctx
        any_replaced = False
        for attr in command_child_attrs(node):
            child = getattr(node, attr)
            if child is not None:
                new_child, replaced = self.walk_close(child, ctx)
                setattr(node, attr, new_child)
                any_replaced = any_replaced or replaced
        # CaseNode: walk cbody in each case
        if isinstance(node, CaseNode):
            for case in node.cases:
                if case.get("cbody") is not None:
                    new_body, replaced = self.walk_close(case["cbody"], ctx)
                    case["cbody"] = new_body
                    any_replaced = any_replaced or replaced
        return NodeResult(ast=node, something_replaced=any_replaced)

    # === Leaf replacement nodes ===

    def visit_pipe(self, node: PipeNode) -> NodeResult:
        return NodeResult(
            ast=node,
            replace_whole=True,
            non_maximal=node.is_background,
            something_replaced=True,
        )

    def visit_command(self, node: CommandNode) -> NodeResult:
        if len(node.arguments) == 0:
            return NodeResult(ast=node, something_replaced=False)
        return NodeResult(ast=node, replace_whole=True, something_replaced=True)

    def visit_background(self, node: BackgroundNode) -> NodeResult:
        return NodeResult(
            ast=node,
            replace_whole=True,
            non_maximal=True,
            something_replaced=True,
        )

    # === No-op nodes (skip children) ===

    def visit_defun(self, node: DefunNode) -> NodeResult:
        return NodeResult(ast=node, something_replaced=False)

    def visit_arith(self, node: ArithNode) -> NodeResult:
        return NodeResult(ast=node, something_replaced=False)
