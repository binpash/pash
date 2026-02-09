"""
Transformation options and state classes for AST preprocessing.
"""

from abc import ABC, abstractmethod
from enum import Enum
import os
import pickle

from shasta.ast_node import AstNode
from shasta.json_to_ast import to_ast_node

from shell_ast.ast_util import string_to_argument, make_command
from parse import from_ast_objects_to_shell
from util import ptempfile


# Runtime executable path - constructed from PASH_TOP environment variable
PASH_TOP = os.environ.get("PASH_TOP", "")
RUNTIME_EXECUTABLE = os.path.join(PASH_TOP, "jit_runtime/jit.sh")


class TransformationType(Enum):
    """Types of AST-to-AST transformations."""

    PASH = "pash"
    AIRFLOW = "airflow"


class AbstractTransformationState(ABC):
    """Base class for transformation state."""

    def __init__(self):
        self._node_counter = 0

    # Node id related
    def get_next_id(self):
        new_id = self._node_counter
        self._node_counter += 1
        return new_id

    def get_current_id(self):
        return self._node_counter - 1

    def get_number_of_ids(self):
        return self._node_counter

    @abstractmethod
    def replace_df_region(
        self, asts, disable_parallel_pipelines=False, ast_text=None
    ) -> AstNode:
        pass


class TransformationState(AbstractTransformationState):
    """Standard PaSh transformation state."""

    def replace_df_region(
        self, asts, disable_parallel_pipelines=False, ast_text=None
    ) -> AstNode:
        ir_filename = ptempfile()

        # Serialize the node in a file
        with open(ir_filename, "wb") as ir_file:
            pickle.dump(asts, ir_file)

        # Serialize the candidate df_region asts back to shell
        # so that the sequential script can be run in parallel to the compilation.
        sequential_script_file_name = ptempfile()
        text_to_output = get_shell_from_ast(asts, ast_text=ast_text)
        with open(sequential_script_file_name, "w", encoding="utf-8") as script_file:
            script_file.write(text_to_output)
        replaced_node = TransformationState.make_call_to_pash_runtime(
            ir_filename, sequential_script_file_name, disable_parallel_pipelines
        )

        return to_ast_node(replaced_node)

    @staticmethod
    def make_call_to_pash_runtime(
        ir_filename, sequential_script_file_name, disable_parallel_pipelines
    ) -> AstNode:
        """
        Make a command that calls the pash runtime with the IR file.
        """
        if disable_parallel_pipelines:
            assignments = [["pash_disable_parallel_pipelines", string_to_argument("1")]]
        else:
            assignments = [["pash_disable_parallel_pipelines", string_to_argument("0")]]
        assignments.append(
            [
                "pash_sequential_script_file",
                string_to_argument(sequential_script_file_name),
            ]
        )
        assignments.append(["pash_input_ir_file", string_to_argument(ir_filename)])

        # Call the runtime
        arguments = [
            string_to_argument("source"),
            string_to_argument(RUNTIME_EXECUTABLE),
        ]
        runtime_node = make_command(arguments, assignments=assignments)
        return runtime_node


class AirflowTransformationState(TransformationState):
    """Airflow transformation state (same as standard PaSh for now)."""

    pass


def get_shell_from_ast(asts, ast_text=None) -> str:
    """Get shell text from AST, using original text if available."""
    if ast_text is None:
        text_to_output = from_ast_objects_to_shell(asts)
    else:
        text_to_output = ast_text
    return text_to_output
