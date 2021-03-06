from definitions.ir.file_id import *

from definitions.ir.dfg_node import *
from ir_utils import string_to_argument, make_command

import config

import os


class RSplit(DFGNode):
    def __init__(self, inputs, outputs, com_name, com_category, com_options = [], 
                 com_redirs = [], com_assignments=[]):
        super().__init__(inputs, outputs, com_name, com_category, 
                         com_options=com_options, 
                         com_redirs=com_redirs, 
                         com_assignments=com_assignments)

## TODO: Make a proper splitter subclass of Node
def make_r_split(input_id, out_ids, r_split_batch_size):
    r_split_bin = os.path.join(config.PASH_TOP, config.config['runtime']['r_split_binary'])
    com_name = Arg(string_to_argument(r_split_bin))
    com_category = "pure"
    com_option = (1, Arg(string_to_argument(str(r_split_batch_size))))
    return RSplit([input_id],
                 out_ids,
                 com_name,
                 com_category,
                 com_options=[com_option])