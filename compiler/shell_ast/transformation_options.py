from abc import ABC, abstractmethod
from enum import Enum
import pickle

from shell_ast.ast_util import *
from shasta.json_to_ast import to_ast_node
from speculative import util_spec
from parse import from_ast_objects_to_shell


## There are two types of ast_to_ast transformations
class TransformationType(Enum):
    PASH = "pash"
    SPECULATIVE = "spec"
    AIRFLOW = "airflow"


class AbstractTransformationState(ABC):
    def __init__(self):
        self._node_counter = 0
        self._loop_counter = 0
        self._loop_contexts = []

    def get_mode(self):
        return TransformationType.PASH

    ## Node id related
    def get_next_id(self):
        new_id = self._node_counter
        self._node_counter += 1
        return new_id

    def get_current_id(self):
        return self._node_counter - 1

    def get_number_of_ids(self):
        return self._node_counter

    ## Loop id related
    def get_next_loop_id(self):
        new_id = self._loop_counter
        self._loop_counter += 1
        return new_id

    def get_current_loop_context(self):
        ## We want to copy that
        return self._loop_contexts[:]

    def get_current_loop_id(self):
        if len(self._loop_contexts) == 0:
            return None
        else:
            return self._loop_contexts[0]

    def enter_loop(self):
        new_loop_id = self.get_next_loop_id()
        self._loop_contexts.insert(0, new_loop_id)
        return new_loop_id

    def exit_loop(self):
        self._loop_contexts.pop(0)

    @abstractmethod
    def replace_df_region(
        self, asts, disable_parallel_pipelines=False, ast_text=None
    ) -> AstNode:
        pass


## Use this object to pass options inside the preprocessing
## trasnformation.
class TransformationState(AbstractTransformationState):
    def replace_df_region(
        self, asts, disable_parallel_pipelines=False, ast_text=None
    ) -> AstNode:
        ir_filename = ptempfile()

        ## Serialize the node in a file
        with open(ir_filename, "wb") as ir_file:
            pickle.dump(asts, ir_file)

        ## Serialize the candidate df_region asts back to shell
        ## so that the sequential script can be run in parallel to the compilation.
        sequential_script_file_name = ptempfile()
        text_to_output = get_shell_from_ast(asts, ast_text=ast_text)
        ## However, if we have the original ast text, then we can simply output that.
        with open(sequential_script_file_name, "w") as script_file:
            script_file.write(text_to_output)
        replaced_node = TransformationState.make_call_to_pash_runtime(
            ir_filename, sequential_script_file_name, disable_parallel_pipelines
        )

        return to_ast_node(replaced_node)

    ## This function makes a command that calls the pash runtime
    ## together with the name of the file containing an IR. Then the
    ## pash runtime should read from this file and continue
    ## execution.
    ##
    ## TODO: At the moment this is written in python but it is in essense a simple shell script.
    ##       Is it possible to make it be a simple string instead of manually creating the AST?
    ##
    ## (MAYBE) TODO: The way I did it, is by calling the parser once, and seeing
    ## what it returns. Maybe it would make sense to call the parser on
    ## the fly to have a cleaner implementation here?
    @staticmethod
    def make_call_to_pash_runtime(
        ir_filename, sequential_script_file_name, disable_parallel_pipelines
    ) -> AstNode:
        ## Disable parallel pipelines if we are in the last command of the script.
        ## ```
        ## pash_disable_parallel_pipelines=1
        ## ```
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

        ## Call the runtime
        arguments = [
            string_to_argument("source"),
            string_to_argument(config.RUNTIME_EXECUTABLE),
        ]
        runtime_node = make_command(arguments, assignments=assignments)
        return runtime_node


## TODO: Turn it into a Transformation State class, and make a subclass for
##       each of the two transformations. It is important for it to be state, because
##       it will need to be passed around while traversing the tree.
class SpeculativeTransformationState(AbstractTransformationState):
    def __init__(self, po_file: str):
        self.partial_order_file = po_file
        self.partial_order_edges = []
        self.partial_order_node_loop_contexts = {}

    def replace_df_region(
        self, asts, disable_parallel_pipelines=False, ast_text=None
    ) -> AstNode:
        text_to_output = get_shell_from_ast(asts, ast_text=ast_text)
        ## Generate an ID
        df_region_id = self.get_next_id()

        ## Get the current loop id and save it so that the runtime knows
        ## which loop it is in.
        loop_id = self.get_current_loop_id()

        ## Determine its predecessors
        ## TODO: To make this properly work, we should keep some state
        ##       in the AST traversal to be able to determine predecessors.
        if df_region_id == 0:
            predecessors = []
        else:
            predecessors = [df_region_id - 1]
        ## Write to a file indexed by its ID
        util_spec.save_df_region(text_to_output, self, df_region_id, predecessors)
        ## TODO: Add an entry point to spec through normal PaSh
        replaced_node = SpeculativeTransformationState.make_call_to_spec_runtime(
            df_region_id, loop_id
        )
        return to_ast_node(replaced_node)

    def get_partial_order_file(self):
        return self.partial_order_file

    def add_edge(self, from_id: int, to_id: int):
        self.partial_order_edges.append((from_id, to_id))

    def get_all_edges(self):
        return self.partial_order_edges

    def add_node_loop_context(self, node_id: int, loop_contexts):
        self.partial_order_node_loop_contexts[node_id] = loop_contexts

    def get_all_loop_contexts(self):
        return self.partial_order_node_loop_contexts

    ## TODO: Make that an actual call to the spec runtime
    @staticmethod
    def make_call_to_spec_runtime(command_id: int, loop_id) -> AstNode:
        assignments = [["pash_spec_command_id", string_to_argument(str(command_id))]]
        if loop_id is None:
            loop_id_str = ""
        else:
            loop_id_str = str(loop_id)

        assignments.append(["pash_spec_loop_id", string_to_argument(loop_id_str)])

        ## Call the runtime
        arguments = [
            string_to_argument("source"),
            string_to_argument(config.RUNTIME_EXECUTABLE),
        ]
        ## Pass all relevant argument to the planner
        runtime_node = make_command(arguments, assignments=assignments)

        return runtime_node


class AirflowTransformationState(TransformationState):
    pass


def get_shell_from_ast(asts, ast_text=None) -> str:
    ## If we don't have the original ast text, we need to unparse the ast
    if ast_text is None:
        text_to_output = from_ast_objects_to_shell(asts)
    else:
        text_to_output = ast_text
    return text_to_output
