from definitions.ir.dfg_node import *

import config

class BigramGMap(DFGNode):
    num_outputs = config.bigram_g_map_num_outputs
    def __init__(self, input, outputs):
        assert(self.num_outputs == 3)
        com_name = Arg(string_to_argument("bigram_aux_map"))
        category = "pure"
        super().__init__([input],
                         outputs, 
                         com_name, 
                         category)
