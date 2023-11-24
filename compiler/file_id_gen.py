from pash_annotations.datatypes.BasicDatatypesWithIOVar import OptionWithIOVar
import config
import os
from definitions.ir.file_id import *
from definitions.ir.nodes.cat import *
from annotations_utils.util_file_descriptors import resource_from_file_descriptor
from pash_annotations.datatypes.BasicDatatypes import ArgStringType
from pash_annotations.datatypes.BasicDatatypesWithIO import (
    FileNameWithIOInfo,
    StdDescriptorWithIOInfo,
    OptionWithIO,
)
from pash_annotations.annotation_generation.datatypes.ParallelizabilityInfo import (
    ParallelizabilityInfo,
)
from pash_annotations.annotation_generation.datatypes.CommandProperties import (
    CommandProperties,
)
from pash_annotations.datatypes.CommandInvocationWithIOVars import (
    CommandInvocationWithIOVars,
)
from annotations_utils.util_parsing import parse_arg_list_to_command_invocation
from annotations_utils.util_cmd_invocations import (
    get_input_output_info_from_cmd_invocation_util,
    get_parallelizability_info_from_cmd_invocation_util,
)
from ir import IR


class FileIdGenerator:
    """
    Generate File IDs
    """

    def __init__(self, next=0, prefix=""):
        assert config.PASH_TMP_PREFIX is not None
        directory = f"{str(uuid.uuid4().hex)}"
        directory_path = os.path.join(config.PASH_TMP_PREFIX, prefix)
        os.makedirs(directory_path)

        self.next = next + 1
        self.prefix = f"{directory}/{prefix}"
        self.dfg_edges = {}
        self.access_map = {}

    def next_file_id(self):
        file_id = FileId(self.next, self.prefix)
        self.next += 1
        return file_id

    def next_temporary_file_id(self):
        file_id = self.next_file_id()
        file_id.make_temporary_file()
        return file_id

    def next_ephemeral_file_id(self):
        file_id = self.next_file_id()
        file_id.make_ephemeral()
        return file_id

    def _create_file_id_for_resource(self, resource):
        """
        Create a file id for a given resource
        """
        file_id = self.next_file_id()
        file_id.set_resource(resource)
        return file_id

    def _add_file_id_vars(self, cmd_with_io):
        new_flagoption_list = self._make_flag_opt_list(cmd_with_io.flag_option_list)

        new_operand_list = self._make_operand_list(cmd_with_io.operand_list)

        if use_stream_in := cmd_with_io.implicit_use_of_streaming_input:
            implicit_stream_in = self._add_var_for_descriptor(use_stream_in)
        else:
            implicit_stream_in = None

        if use_stream_out := cmd_with_io.implicit_use_of_streaming_output:
            implicit_stream_out = self._add_var_for_descriptor(use_stream_out)
        else:
            implicit_stream_out = None

        command_invocation_with_io_vars = CommandInvocationWithIOVars(
            cmd_name=cmd_with_io.cmd_name,
            flag_option_list=new_flagoption_list,
            operand_list=new_operand_list,
            implicit_use_of_streaming_input=implicit_stream_in,
            implicit_use_of_streaming_output=implicit_stream_out,
            access_map=self.access_map,
        )
        return command_invocation_with_io_vars, self.dfg_edges

    def _make_operand_list(self, operand_list):
        new_op_list = []
        for operand in operand_list:
            if isinstance(operand, FileNameWithIOInfo) or isinstance(
                operand, StdDescriptorWithIOInfo
            ):
                fid_id = self._add_var_for_descriptor(operand)
                new_op_list.append(fid_id)
            else:
                new_op_list.append(operand)
        return new_op_list

    def _make_flag_opt_list(self, opt_list):
        new_opt_list = []
        for flag_opt in opt_list:
            if isinstance(flag_opt, OptionWithIO) and not isinstance(
                flag_opt.option_arg, ArgStringType
            ):
                fid_id = self._add_var_for_descriptor(flag_opt.option_arg)
                new_option = OptionWithIOVar(flag_opt.name, fid_id)
                new_opt_list.append(new_option)
            else:  # Flag
                new_opt_list.append(flag_opt)
        return new_opt_list

    def _add_var_for_descriptor(self, operand):
        resource = resource_from_file_descriptor(operand)
        file_id = self._create_file_id_for_resource(resource)
        fid_id = file_id.get_ident()
        self.dfg_edges[fid_id] = (file_id, None, None)
        self.access_map[fid_id] = operand.get_access()
        return fid_id

    def compile_command_to_DFG(self, command, options, redirections=[]):
        command_invocation = parse_arg_list_to_command_invocation(command, options)
        io_info = get_input_output_info_from_cmd_invocation_util(command_invocation)
        if io_info is None:
            raise Exception(
                f"InputOutputInformation for {format_arg_chars(command)} not provided so considered side-effectful."
            )
        elif io_info.has_other_outputs():
            raise Exception(
                f"Command {format_arg_chars(command)} has outputs other than streaming."
            )
        para_info = get_parallelizability_info_from_cmd_invocation_util(
            command_invocation
        )
        if para_info is None:
            # defaults to no parallelizer's and all properties False
            para_info = ParallelizabilityInfo()
        cmd_with_io = io_info.apply_input_output_info_to_command_invocation(
            command_invocation
        )

        parallelizer_lst, round_robin, commutative = para_info.unpack_info()
        property_dict = [
            {
                "round_robin_compatible_with_cat": round_robin,
                "is_commutative": commutative,
            }
        ]
        cmd_related_properties = CommandProperties(property_dict)

        ## TODO: Make an empty IR and add edges and nodes incrementally (using the methods defined in IR).

        ## Add all inputs and outputs to the DFG edges
        cmd_with_io_vars, dfg_edges = self._add_file_id_vars(cmd_with_io)
        com_redirs = redirections
        ## TODO: Add assignments
        com_assignments = []

        ## Assume: Everything must be completely expanded
        ## TODO: Add an assertion about that.
        dfg_node = DFGNode(
            cmd_with_io_vars,
            com_redirs=com_redirs,
            com_assignments=com_assignments,
            parallelizer_list=parallelizer_lst,
            cmd_related_properties=cmd_related_properties,
        )
        # log(f'Dfg node: {dfg_node}')
        node_id = dfg_node.get_id()

        ## Assign the from, to node in edges
        for fid_id in dfg_node.get_input_list():
            fid, from_node, to_node = dfg_edges[fid_id]
            assert to_node is None
            dfg_edges[fid_id] = (fid, from_node, node_id)

        for fid_id in dfg_node.get_output_list():
            fid, from_node, to_node = dfg_edges[fid_id]
            assert from_node is None
            dfg_edges[fid_id] = (fid, node_id, to_node)

        dfg_nodes = {node_id: dfg_node}
        dfg = IR(dfg_nodes, dfg_edges)
        return dfg
