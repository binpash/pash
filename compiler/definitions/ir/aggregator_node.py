from definitions.ir.dfg_node import *

## This class corresponds to a generic n-ary aggregator
##
## TODO: Do we need to do anything special for binary aggregators?
class AggregatorNode(DFGNode):
    def __init__(self, old_node, input_ids, output_id):

        ## Check if an aggregator can be instantiated from the node
        if(old_node.com_aggregator is None):
            log("Error: Node:", old_node, "does not contain information to instantiate an aggregator!")
            exit(1)

        ## The name of the aggregator command
        agg_name_string = old_node.com_aggregator.name
        name = Arg(string_to_argument(agg_name_string))

        ## TODO: The category should also be acquired through annotations (and maybe should be asserted to be at most pure)
        com_category="pure"

        ## TODO: Not sure if redirections need to be copied to new function.
        com_redirs = [redir.to_ast() for redir in old_node.com_redirs]
        super().__init__(input_ids,
                         [output_id], 
                         name,
                         com_category, 
                         com_options=old_node.com_options, 
                         com_redirs=com_redirs, 
                         com_assignments=old_node.com_assignments)
        
        ## TODO: This assumes that all options from the old function are copied to the new.
        ##
        ## TODO: If we need a behavior where we don't keep the old flags, we can extend this
        self.append_options(old_node.com_aggregator.options)

        log("Generic Aggregator Created:", self)
