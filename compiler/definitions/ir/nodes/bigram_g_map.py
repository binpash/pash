from definitions.ir.dfg_node import *

class BigramGMap(DFGNode):
    def __init__(self, input, outputs):
        com_name = Arg(string_to_argument("bigram_aux_map"))
        category = "pure"
        super().__init__([input],
                         outputs, 
                         com_name, 
                         category)
