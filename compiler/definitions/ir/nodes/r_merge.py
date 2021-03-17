from definitions.ir.dfg_node import *

class RMerge(DFGNode):
    def __init__(self, inputs, outputs, com_name, com_category,
                 com_options = [], com_redirs = [], com_assignments=[]):
        super().__init__(inputs, outputs, com_name, com_category,
                         com_options=com_options, 
                         com_redirs=com_redirs, 
                         com_assignments=com_assignments)

def make_r_merge_node(inputs, output):
    r_merge_bin = os.path.join(config.PASH_TOP, config.config['runtime']['r_merge_binary'])
    com_name = Arg(string_to_argument(r_merge_bin))
    com_category = "pure"
    return RMerge(inputs,
                  [output],
                  com_name, 
                  com_category)
