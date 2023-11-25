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
from pash_annotations.datatypes.CommandInvocationWithIOVars import (
    CommandInvocationWithIOVars,
)


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

    def add_file_id_vars(self, cmd_with_io):
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

    def _create_file_id_for_resource(self, resource):
        """
        Create a file id for a given resource
        """
        file_id = self.next_file_id()
        file_id.set_resource(resource)
        return file_id

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
