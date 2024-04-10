import os
from definitions.ir.dfg_node import *


class DFSSplitReader(DFGNode):
    def __init__(self, inputs, outputs, com_name, com_category,
                 com_options=[], com_redirs=[], com_assignments=[]):

        super().__init__(inputs, outputs, com_name, com_category,
                         com_options=com_options,
                         com_redirs=com_redirs,
                         com_assignments=com_assignments)

    def set_server_address(self, addr):  # ex addr: 127.0.0.1:50051
        self.com_options.append((3, Arg(string_to_argument(f"--addr {addr}"))))


def make_dfs_split_reader_node(inputs, output, split_num, subblock_num, subblock_cnt):
    split_reader_bin = os.path.join(
        config.DISH_TOP, config.config['runtime']['dfs_split_reader_binary'])
    com_name = Arg(string_to_argument(split_reader_bin))
    com_category = "pure"
    options = []
    options.append((1, Arg(string_to_argument(f"--split {split_num}"))))
    options.append((2, Arg(string_to_argument(f"--subSplit {subblock_num}"))))
    options.append((3, Arg(string_to_argument(f"--numSplits {subblock_cnt}"))))

    return DFSSplitReader(inputs,
                          [output],
                          com_name,
                          com_category,
                          options)
