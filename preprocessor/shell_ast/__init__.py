"""
Shell AST manipulation for preprocessing.
"""

from shell_ast.ast_to_ast import preprocess_node, replace_ast_regions
from shell_ast.transformation_options import (
    TransformationType,
    TransformationState,
    SpeculativeTransformationState,
    AirflowTransformationState,
)

__all__ = [
    "preprocess_node",
    "replace_ast_regions",
    "TransformationType",
    "TransformationState",
    "SpeculativeTransformationState",
    "AirflowTransformationState",
]
