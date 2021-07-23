from definitions.ir.dfg_node import *

## This class corresponds to a generic n-ary aggregator
##
## TODO: Do we need to do anything special for binary aggregators?
class MapperAggregatorNode(DFGNode):
    def __init__(self, old_node, input_ids, output_ids, name_string, new_options):

        ## The name of the aggregator command
        name = Arg(string_to_argument(name_string))

        ## TODO: The category should also be acquired through annotations (and maybe should be asserted to be at most pure)
        com_category="pure"

        ## TODO: Not sure if redirections need to be copied to new function.
        com_redirs = [redir.to_ast() for redir in old_node.com_redirs]
        super().__init__(input_ids,
                         output_ids, 
                         name,
                         com_category, 
                         com_options=old_node.com_options, 
                         com_redirs=com_redirs, 
                         com_assignments=old_node.com_assignments)
        
        ## TODO: This assumes that all options from the old function are copied to the new.
        ##
        ## TODO: If we need a behavior where we don't keep the old flags, we can extend this
        self.append_options(new_options)


class AggregatorNode(MapperAggregatorNode):
    def __init__(self, old_node, input_ids, output_ids):

        ## Check if an aggregator can be instantiated from the node
        if(old_node.com_aggregator is None):
            log("Error: Node:", old_node, "does not contain information to instantiate an aggregator!")
            exit(1)

        ## The name of the aggregator command
        agg_name_string = old_node.com_aggregator.name
        new_options = old_node.com_aggregator.options

        super().__init__(old_node, input_ids, output_ids, agg_name_string, new_options)

        log("Generic Aggregator Created:", self)

class MapperNode(MapperAggregatorNode):
    def __init__(self, old_node, input_ids, output_ids):

        ## Check if an mapper can be instantiated from the node
        if(old_node.com_mapper is None):
            log("Error: Node:", old_node, "does not contain information to instantiate a mapper!")
            exit(1)

        ## The name of the aggregator command
        mapper_name_string = old_node.com_mapper.name
        new_options = old_node.com_mapper.options

        super().__init__(old_node, input_ids, output_ids, mapper_name_string, new_options)

        log("Generic Mapper Created:", self)
