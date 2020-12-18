from command_categories import *
from util import *
from ir_utils import *

from definitions.ir.redirection import *
from definitions.ir.resource import *

## Assumption: Everything related to a DFGNode must be already expanded.
## TODO: Ensure that this is true with assertions
class DFGNode:
    ## TODO: Make inputs + outputs be structure of fids instead of list. 
    ##       This will allow us to handle comm and other commands with static inputs. 
    ##
    ## inputs : list of fid_ids (that can be used to retrieve fid from edges)
    ## outputs : list of fid_ids 
    ## com_name : command name Arg
    ## com_options : list of tuples with the option index and the argument Arg
    ## com_redirs : list of redirections
    ## com_assignments : list of assignments
    def __init__(self, inputs, outputs, com_name, com_category, com_options = [], 
                 com_redirs = [], com_assignments=[]):
        self.inputs = inputs
        self.outputs = outputs
        self.com_name = com_name
        self.com_category = com_category
        self.com_options = com_options
        self.com_redirs = [Redirection(redirection) for redirection in com_redirs]
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

    def get_input_fids(self, edges):
        return [fid for _, fid in self.get_input_ids_fids()]

    def get_output_fids(self, edges):
        return [fid for _, fid in self.get_output_ids_fids()]

    def get_input_ids_fids(self, edges):
        return [(input_edge_id, edges[input_edge_id][0]) for input_edge_id in self.inputs]

    def get_output_ids_fids(self, edges):
        return [(output_edge_id, edges[output_edge_id][0]) for output_edge_id in self.outputs]

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
            ##
            ## TODO: If one of the redirected outputs or inputs is changed in the IR 
            ##       (e.g. `cat < s1` was changed to read from an ephemeral file `cat < "#file5"`)
            ##       this needs to be changed in the redirections too. Maybe we can modify redirections
            ##       when replacing fid.
            ##
            ## It seems that if we have already applied redirections we might not need to
            ## care about them anymore (since they will be created in new_redirs.)
            ##
            ## redirs = [redir.to_ast() for redir in self.com_redirs]
            ##
            ## At the moment we do not reprint redirections here (we only produce redirections
            ## where we recreate arguments and redirections).
            redirs = []
            assignments = self.com_assignments
            ## Start filling in the arguments
            opt_arguments = []
            for i, opt in self.com_options:
                ## Pad the argument list with None 
                opt_arguments = pad(opt_arguments, i)
                opt_arguments[i] = opt.to_ast()

            com_name_ast = self.com_name.to_ast()
            option_asts = [opt.to_ast() for _, opt in self.com_options]

            ##
            ## 1. Find the input and output fids
            ## 2. Construct the rest of the arguments and input/output redirections according to
            ##    the command IO
            input_fids = [edges[in_id][0] for in_id in self.inputs]
            output_fids = [edges[out_id][0] for out_id in self.outputs]
            rest_argument_fids, new_redirs = create_command_arguments_redirs(com_name_ast,
                                                                             option_asts,
                                                                             input_fids,
                                                                             output_fids)
            
            ## Transform the rest of the argument fids to arguments
            rest_arguments = [fid.to_ast() for fid in rest_argument_fids]

            ## Interleave the arguments since options args might contain gaps.
            arguments = interleave_args(opt_arguments, rest_arguments) 

            all_arguments = [com_name_ast] + arguments
            all_redirs = redirs + new_redirs

            node = make_command(all_arguments, redirections=all_redirs, assignments=assignments)
        return node

    ## This method applies the redirections to get the correct, inputs, outputs of a node.
    ##
    ## WARNING: For now it only works with 'To' redirections for
    ## stdout, and it applies them by adding a resource to the stdout
    ## of the command. It also keeps them for possible future usage.
    ##
    ## TODO: Properly handle all redirections. This requires a nice
    ## abstraction. Maybe the best way would be to keep them around
    ## and always recompute inputs/outputs when needed by following
    ## the redirections.
    ##
    ## TODO: Is it correct to apply redirections from left to right?
    def apply_redirections(self, edges):
        unhandled_redirs = []
        for redirection in self.com_redirs:
            ## Handle To redirections that have to do with stdout
            if (redirection.is_to_file() and redirection.is_for_stdout()):
                # log(redirection)
                file_resource = FileResource(redirection.file_arg)
                success = False
                for i in range(len(self.outputs)):
                    output_edge_id = self.outputs[i]
                    output_fid = edges[output_edge_id][0]
                    if(output_fid.has_file_descriptor_resource()
                       and output_fid.resource.is_stdout()):
                        success = True
                        edges[output_edge_id][0].set_resource(file_resource)
                        # self.outputs[i].set_resource(file_resource)
                assert(success)
            elif (redirection.is_from_file() and redirection.is_for_stdin()):
                # log(redirection)
                file_resource = FileResource(redirection.file_arg)
                success = False
                for i in range(len(self.inputs)):
                    input_edge_id = self.inputs[i]
                    input_fid = edges[input_edge_id][0]
                    if(input_fid.has_file_descriptor_resource()
                       and input_fid.resource.is_stdin()):
                        success = True
                        edges[input_edge_id][0].set_resource(file_resource)
                assert(success)
            else:
                log("Warning -- Unhandled redirection:", redirection)
                unhandled_redirs.append(redirection)
                ## TODO: I am not sure if this is the correct way to handle unhandled redirections.
                ##       Does it make any sense to keep them and have them in the Final AST.
                raise NotImplementedError()


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