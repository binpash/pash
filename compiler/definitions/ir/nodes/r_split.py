import os

import config

from definitions.ir.dfg_node import *
from definitions.ir.file_id import *
from ir_utils import string_to_argument, make_command

class RSplit(DFGNode):
    def __init__(self, inputs, outputs, com_name, com_category, com_options = [], 
                 com_redirs = [], com_assignments=[]):
        super().__init__(inputs, outputs, com_name, com_category, 
                         com_options=com_options, 
                         com_redirs=com_redirs, 
                         com_assignments=com_assignments)
    
    ## TODO: Generalize this code (for this and SortGReduce) to be able to add an option to any command.
    def add_r_flag(self):
        assert(len(self.com_options) <= 1)
            
        ## Add -r in r_split
        new_opt = (0, Arg(string_to_argument("-r")))
        shifted_options = [(i+1, opt) for i, opt in self.com_options]
        self.com_options = [new_opt] + shifted_options

    ## This is not a proper option check. It just works if the r_flag is added as a separate option.
    def has_r_flag(self):
        option_strings = [str(opt) for i, opt in self.com_options]
        return ("-r" in option_strings)


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