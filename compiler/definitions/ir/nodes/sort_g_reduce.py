from definitions.ir.dfg_node import *

class SortGReduce(DFGNode):
    def __init__(self, old_node, ids):
        assert(str(old_node.com_name) == "sort")
        input_ids = ids[:-1]
        output_id = ids[-1]
        com_category="pure"
        if(len(old_node.com_options) > 0):
            max_opt_index = max([i for i, _opt in old_node.com_options])
        else:
            max_opt_index = -1
        com_options = old_node.com_options + [(max_opt_index+1, Arg(string_to_argument("-m")))]
        com_redirs = [redir.to_ast() for redir in old_node.com_redirs]
        super().__init__(input_ids,
                         [output_id], 
                         old_node.com_name, 
                         com_category, 
                         com_options=com_options, 
                         com_redirs=com_redirs, 
                         com_assignments=old_node.com_assignments)
