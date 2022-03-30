from definitions.ir.dfg_node import *

class HDFSCat(DFGNode):
    def __init__(self, inputs, outputs, com_name, com_category,
                 com_options = [], com_redirs = [], com_assignments=[]):
        assert(str(com_name) == "hdfs")
        assert(str(com_options[0][1]) == "dfs" and str(com_options[1][1]) == "-cat") 
        super().__init__(inputs, outputs, com_name, com_category,
                         com_options=com_options, 
                         com_redirs=com_redirs, 
                         com_assignments=com_assignments)
