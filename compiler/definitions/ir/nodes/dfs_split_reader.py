import os
from definitions.ir.dfg_node import *
from dspash.hdfs_utils import HDFSFileConfig


class DFSSplitReader(DFGNode):
    def __init__(self, inputs, outputs, com_name, com_category,
                 com_options=[], com_redirs=[], com_assignments=[]):

        super().__init__(inputs, outputs, com_name, com_category,
                         com_options=com_options,
                         com_redirs=com_redirs,
                         com_assignments=com_assignments)


def make_dfs_split_reader_node(inputs, output, split_num, subblock_num, subblock_cnt, ft, file_config: HDFSFileConfig):
    split_reader_bin = os.path.join(config.DISH_TOP, config.config['runtime']['dfs_split_reader_binary'])
    com_name = Arg(string_to_argument(split_reader_bin))
    com_category = "pure"
    options = []
    options.append((1, Arg(string_to_argument(f"--split {split_num}"))))
    options.append((2, Arg(string_to_argument(f"--ft {ft}"))))
    if ft == "optimized" or ft == "dynamic":
        options.append((3, Arg(string_to_argument(f"--subblockNum {subblock_num}"))))
        options.append((4, Arg(string_to_argument(f"--subblockCnt {subblock_cnt}"))))
        options.append((5, Arg(string_to_argument(f"--path=\\''{file_config.blocks[split_num].path}'\\'"))))
        if split_num < len(file_config.blocks) - 1 and subblock_num == subblock_cnt - 1:
            options.append((6, Arg(string_to_argument(f"--nextBlockPath=\\''{file_config.blocks[split_num + 1].path}'\\'"))))
            options.append((7, Arg(string_to_argument(f"--nextBlockHosts={','.join(file_config.blocks[split_num + 1].hosts)}"))))

    return DFSSplitReader(inputs,
                          [output],
                          com_name,
                          com_category,
                          options)
