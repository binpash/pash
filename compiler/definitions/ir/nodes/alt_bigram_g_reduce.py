from definitions.ir.dfg_node import *

class AltBigramGReduce(DFGNode):
    def __init__(self, old_node, edge_ids):
        assert(str(old_node.com_name) == "alt_bigrams_aux")
        com_name = Arg(string_to_argument("alt_bigram_aux_reduce"))
        com_category = "pure"
        super().__init__(edge_ids[:2],
                         edge_ids[2:],
                         com_name,
                         com_category)
