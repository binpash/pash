from command_categories import *
from util import *
from ir_utils import *

## Assumption: Everything related to a DFGNode must be already expanded.
## TODO: Ensure that this is true with assertions
class DFGNode:
    ## TODO: Make inputs + outputs be structure of fids instead of list. 
    ##       This will allow us to handle comm and other commands with static inputs. 
    ##
    ## inputs : list of fids
    ## outputs : list of fids
    ## com_name : command name Arg
    ## com_options : list of arguments Arg
    ## com_redirs : list of redirections
    ## com_assignments : list of assignments
    def __init__(self, inputs, outputs, com_name, com_category, com_options = [], 
                 com_redirs = [], com_assignments=[]):
        self.inputs = inputs
        self.outputs = outputs
        self.com_name = com_name
        self.com_category = com_category
        self.com_options = com_options
        self.com_redirs = com_redirs
        self.com_assignments = com_assignments    

    def __repr__(self):
        prefix = "Node"
        if (self.com_category == "stateless"):
            prefix = "Stateless"
        elif (self.com_category == "pure"):
            prefix = "Pure"
        output = "{}: \"{}\" in:{} out:{}".format(
            prefix, self.com_name, 
            self.inputs,
            self.outputs)
        return output

    ## TODO: Replace DFGNodes in the nodes of the IR
    ##       - During compilation (AST to IR) create DFGNodes
    ##       - Modify all IR methods to work with DFGNodes
    ##       - Remove UnionFind from files now that we don't have to have references to arguments.
    ##       - Remove flatten/unflatten for Fids
    ##       - Write functions that are able to recreate ASTNode from DFGNode

    def is_at_most_pure(self):
        return (self.com_category in ["stateless", "pure", "parallelizable_pure"])


    def is_pure_parallelizable(self):
        return (self.com_category == "parallelizable_pure" or
                (self.com_category == "pure"
                 and str(self.com_name) in list(map(get_command_from_definition,
                                                    config.parallelizable_pure_commands))))

    ## TODO: Improve this functio to be separately implemented for different special nodes,
    ##       such as cat, eager, split, etc...
    def to_ast(self, edges, drain_streams):    
        ## TODO: We might not want to implement this at all actually
        if (drain_streams):
            raise NotImplementedError()
        else:
            ## TODO: Properly handle redirections
            redirs = self.com_redirs
            assignments = self.com_assignments
            option_asts = [opt.to_ast() for opt in self.com_options]
            arguments = [self.com_name.to_ast()] + option_asts 

            ##
            ## 1. Find the input and output fids
            ## 2. Construct the rest of the arguments and input/output redirections according to
            ##    the command IO
            input_fids = [edges[in_id][0] for in_id in self.inputs]
            output_fids = [edges[out_id][0] for out_id in self.outputs]
            rest_argument_fids, new_redirs = create_command_arguments_redirs(self.com_name.to_ast(),
                                                                         option_asts,
                                                                         input_fids,
                                                                         output_fids)
            
            ## Transform the rest of the argument fids to arguments
            rest_arguments = [fid.to_ast() for fid in rest_argument_fids] 

            all_arguments = arguments + rest_arguments
            all_redirs = redirs + new_redirs

            node = make_command(all_arguments, redirections=all_redirs, assignments=assignments)
        return node



    ## This renames the from_id (wherever it exists in inputs ot outputs)
    ## to the to_id.
    ##
    ## TODO: Make sure we don't need to change redirections here.
    def replace_edge(self, from_id, to_id):
        new_inputs = []
        for input_id in self.inputs:
            if(input_id == from_id):
                new_input_id = to_id
            else:
                new_input_id = input_id
            new_inputs.append(new_input_id)

        new_outputs = []
        for output_id in self.outputs:
            if(output_id == from_id):
                new_output_id = to_id
            else:
                new_output_id = output_id
            new_outputs.append(new_output_id)
        self.inputs = new_inputs
        self.outputs = new_outputs