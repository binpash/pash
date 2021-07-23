from definitions.ir.dfg_node import *

## TODO: Make that generic so that aggregators can be identified using annotations.
class SortGReduce(DFGNode):
    def __init__(self, old_node, ids):

        ## TODO: The name of the aggregator command can be extracted from the annotation
        name = Arg(string_to_argument("sort"))

        ## TODO: This assumes that the aggregator is n-ary, but we might only want to support generic binary aggregators for start
        input_ids = ids[:-1]
        output_id = ids[-1]

        ## TODO: The category should also be acquired through annotations (and maybe should be asserted to be at most pure)
        com_category="pure"

        ## TODO: This assumes that all options from the old function are copied to the new.
        if(len(old_node.com_options) > 0):
            max_opt_index = max([i for i, _opt in old_node.com_options])
        else:
            max_opt_index = -1
        com_options = old_node.com_options + [(max_opt_index+1, Arg(string_to_argument("-m")))]

        ## TODO: Not sure if redirections need to be copied to new function.
        com_redirs = [redir.to_ast() for redir in old_node.com_redirs]
        super().__init__(input_ids,
                         [output_id], 
                         name,
                         com_category, 
                         com_options=com_options, 
                         com_redirs=com_redirs, 
                         com_assignments=old_node.com_assignments)
