import copy
from command_categories import *
from util import *
from ir_utils import *

from definitions.ir.redirection import *
from definitions.ir.resource import *

import config

## Assumption: Everything related to a DFGNode must be already expanded.
## TODO: Ensure that this is true with assertions
class DFGNode:
    ## Unique identifier for nodes
    next_id = 0

    ## inputs : tuple of lists of fid_ids (that can be used to retrieve fid from edges)
    ## outputs : list of fid_ids 
    ## com_name : command name Arg
    ## com_category : string denoting category
    ## input_consumption_mode : enumeration
    ## com_properties : properties such as commutativity
    ## com_aggregator : a class that contains necessary information to instantiate an aggregator
    ## com_options : list of tuples with the option index and the argument Arg
    ## com_redirs : list of redirections
    ## com_assignments : list of assignments
    def __init__(self, inputs, outputs, com_name, com_category,
                 com_properties = [],
                 com_aggregator = None,
                 com_options = [],
                 com_redirs = [],
                 com_assignments=[]):
        ## Add a unique identifier to each DFGNode since id() is not guaranteed to be unique for objects that have different lifetimes.
        ## This leads to issues when nodes are deleted and new ones are created, leading to id() clashes between them
        self.id = DFGNode.next_id
        DFGNode.next_id += 1

        self.set_inputs(inputs)
        self.outputs = outputs
        self.com_name = com_name
        self.com_category = com_category
        self.com_properties = com_properties
        self.com_aggregator = com_aggregator
        self.com_options = com_options
        self.com_redirs = [Redirection(redirection) for redirection in com_redirs]
        self.com_assignments = com_assignments

        # log("Node created:", self.id, self)

    def __repr__(self):
        prefix = "Node"
        if (self.com_category == "stateless"):
            prefix = "Stateless"
        elif (self.com_category == "pure"):
            prefix = "Pure"
        if (self.is_commutative()):
            prefix = 'Commutative ' + prefix
        output = "{}: \"{}\" in:{} out:{}".format(
            prefix, self.com_name, 
            self.get_input_list(),
            self.outputs)
        return output

    def get_id(self):
        return self.id

    ## Copying requires setting the id to a new one too
    def copy(self):
        node_copy = copy.deepcopy(self)
        node_copy.id = DFGNode.next_id
        DFGNode.next_id += 1
        return node_copy

    ## TODO: Make that a proper class.
    def set_inputs(self, inputs):
        if(isinstance(inputs, list)):
            self.inputs = ([], inputs)
        elif(isinstance(inputs, tuple)):
            self.inputs = inputs
        else:
            raise NotImplementedError()

    def get_input_list(self):
        return (self.inputs[0] + self.inputs[1])
    
    def get_standard_inputs(self):
        return self.inputs[1]
    
    def get_configuration_inputs(self):
        return self.inputs[0]

    def is_at_most_pure(self):
        return (self.com_category in ["stateless", "pure", "parallelizable_pure"])

    def is_parallelizable(self):
        return (self.is_pure_parallelizable() or self.is_stateless())

    def is_stateless(self):
        return (self.com_category == "stateless")

    def is_pure_parallelizable(self):
        return (self.com_category == "parallelizable_pure" or
                (self.com_category == "pure"
                 and str(self.com_name) in list(map(get_command_from_definition,
                                                    config.parallelizable_pure_commands))))

    def is_commutative(self):
        return ('commutative' in self.com_properties)

    ## kk: 2021-07-23 Not totally sure if that is generally correct. Tests will say ¯\_(ツ)_/¯
    ##     I think it assumes that new options can be added in the beginning if there are no options already
    def append_options(self, new_options):
        if(len(self.com_options) > 0):
            max_opt_index = max([i for i, _opt in self.com_options])
        else:
            max_opt_index = -1
        new_com_options = [(max_opt_index + 1 + i, Arg(string_to_argument(opt))) 
                           for i, opt in enumerate(new_options)]
        self.com_options = self.com_options + new_com_options

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
            input_fids = [edges[in_id][0] for in_id in self.get_input_list()]
            output_fids = [edges[out_id][0] for out_id in self.outputs]
            rest_argument_fids, new_redirs = create_command_arguments_redirs(com_name_ast,
                                                                             option_asts,
                                                                             input_fids,
                                                                             output_fids)
            
            ## Transform the rest of the argument fids to arguments
            ## Since some of the rest_arguments can be None (they only contain inputs and outputs)
            ## we need to make sure that we don't turn None objects to asts.
            ##
            ## The None fields need to be filtered out because they are taken care of by the interleave function.
            ##
            ## TODO: Is this actually OK?
            rest_arguments = [fid.to_ast()
                              for fid in rest_argument_fids
                              if not fid is None]

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
                for input_edge_id in self.get_input_list():
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
    ##
    ## TODO: Make this a method of graph to change the from, to too.
    def replace_edge(self, from_id, to_id):
        new_config_inputs = self.replace_edge_in_list(self.inputs[0], from_id, to_id)
        new_standard_inputs = self.replace_edge_in_list(self.inputs[1], from_id, to_id)
        new_outputs = self.replace_edge_in_list(self.outputs, from_id, to_id)
        
        self.set_inputs((new_config_inputs, new_standard_inputs))
        self.outputs = new_outputs

    ## TODO: There must be a lib function to do this.
    def replace_edge_in_list(self, edge_ids, from_id, to_id):
        new_edge_ids = []
        for id in edge_ids:
            if(id == from_id):
                new_edge_id = to_id
            else:
                new_edge_id = id
            new_edge_ids.append(new_edge_id)
        return new_edge_ids

    ## Get the file names of the outputs of the map commands. This
    ## differs if the command is stateless, pure that can be
    ## written as a map and a reduce, and a pure that can be
    ## written as a generalized map and reduce.
    def get_map_output_files(self, input_edge_ids, fileIdGen):
        assert(self.is_parallelizable())
        if(self.com_category == "stateless"):
            map_output_fids = [fileIdGen.next_ephemeral_file_id() for in_fid in input_edge_ids]
        elif(self.is_pure_parallelizable()):
            map_output_fids = self.pure_get_map_output_files(input_edge_ids, fileIdGen)
        else:
            log("Unreachable code reached :(")
            assert(False)
            ## This should be unreachable
        
        return map_output_fids

    ## TODO: Fix this somewhere in the annotations and not in the code
    def pure_get_map_output_files(self, input_edge_ids, fileIdGen):
        assert(self.is_pure_parallelizable())
        if(str(self.com_name) == "sort"):
            new_output_fids = [[fileIdGen.next_ephemeral_file_id()] for in_fid in input_edge_ids]
        elif(str(self.com_name) == "custom_sort"):
            new_output_fids = [[fileIdGen.next_ephemeral_file_id()] for in_fid in input_edge_ids]
        elif(str(self.com_name) == "bigrams_aux"):
            new_output_fids = [[fileIdGen.next_ephemeral_file_id()
                                for i in range(config.bigram_g_map_num_outputs)]
                               for in_fid in input_edge_ids]
        elif(str(self.com_name) == "alt_bigrams_aux"):
            new_output_fids = [[fileIdGen.next_ephemeral_file_id()] for in_fid in input_edge_ids]
        elif(str(self.com_name) == "uniq"):
            new_output_fids = [[fileIdGen.next_ephemeral_file_id()] for in_fid in input_edge_ids]
        else:
            log("Error: Map outputs for command:", self.com_name, "were not found!")
            raise NotImplementedError()
        return new_output_fids
