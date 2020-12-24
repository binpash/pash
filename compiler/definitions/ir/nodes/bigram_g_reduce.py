from definitions.ir.dfg_node import *

class BigramGReduce(DFGNode):
    def __init__(self, old_node, edge_ids):
        assert(str(old_node.com_name) == "bigrams_aux")
        com_name = Arg(string_to_argument("bigram_aux_reduce"))
        com_category = "pure"
        super().__init__(edge_ids[:6],
                         edge_ids[6:],
                         com_name,
                         com_category)
