import os
from definitions.ir.dfg_node import *


class DFSSplitReader(DFGNode):
    def __init__(
        self,
        inputs,
        outputs,
        com_name,
        com_category,
        com_options=None,
        com_redirs=None,
        com_assignments=None,
    ):
        com_options = [] if com_options is None else com_options
        com_redirs = [] if com_redirs is None else com_redirs
        com_assignments = [] if com_assignments is None else com_assignments

        super().__init__(
            inputs,
            outputs,
            com_name,
            com_category,
            com_options=com_options,
            com_redirs=com_redirs,
            com_assignments=com_assignments,
        )

    def set_server_address(self, addr):  # ex addr: 127.0.0.1:50051
        self.com_options.append((3, Arg.string_to_arg(f"--addr {addr}")))


def make_dfs_split_reader_node(inputs, output, split_num, prefix):
    split_reader_bin = os.path.join(
        config.PASH_TOP, config.config["runtime"]["dfs_split_reader_binary"]
    )
    com_name = Arg.string_to_arg(split_reader_bin)
    com_category = "pure"
    options = []
    options.append((1, Arg.string_to_arg(f"--prefix '{prefix}'")))
    options.append((2, Arg.string_to_arg(f"--split {split_num}")))

    return DFSSplitReader(inputs, [output], com_name, com_category, options)
