import copy
from command_categories import *
from definitions.ir.redirection import *
from definitions.ir.resource import *

from annotations_utils.util_cmd_invocations import to_node_cmd_inv_with_io_vars, construct_property_container_from_list_of_properties

import sys
from config import get_path_annotation_repo
sys.path.insert(1, get_path_annotation_repo())

from util import return_empty_list_if_none_else_itself, return_default_if_none_else_itself

## Assumption: Everything related to a DFGNode must be already expanded.
## TODO: Ensure that this is true with assertions
class DFGNode:
    ## Unique identifier for nodes
    next_id = 0

    ## cmd_invocation_with_io_vars : command invocation data structure with edge ids as symbolic variables for filenames etc.
    ## com_redirs : list of redirections
    ## com_assignments : list of assignments
    ## parallelizer_list : list of parallelizers for this DFGNode
    ## cmd_related_properties : dict to store properties like commutativity
    def __init__(self,
                 cmd_invocation_with_io_vars,
                 com_redirs = [],
                 com_assignments=[],
                 parallelizer_list=None,
                 cmd_related_properties=None,
                 ):
        # TODO []: default parameters!

        ## @KK: can this be deleted? Was there another id in the member attributes before?
        ## Add a unique identifier to each DFGNode since id() is not guaranteed to be unique for objects that have different lifetimes.
        ## This leads to issues when nodes are deleted and new ones are created, leading to id() clashes between them
        self.id = DFGNode.next_id
        DFGNode.next_id += 1

        self.com_redirs = [Redirection(redirection) for redirection in com_redirs]
        self.com_assignments = com_assignments
        self.parallelizer_list = return_empty_list_if_none_else_itself(parallelizer_list)
        default_cmd_properties = construct_property_container_from_list_of_properties([])
        self.cmd_related_properties = return_default_if_none_else_itself(cmd_related_properties, default_cmd_properties)
        self.cmd_invocation_with_io_vars = cmd_invocation_with_io_vars
        # log("Node created:", self.id, self)

    def __repr__(self):
        # TODO: add other attributes
        return str(self.cmd_invocation_with_io_vars)

    ## Generates a dot node for the DFG node
    def add_dot_node(self, dot, node_id):
        label = self.get_dot_label()
        dot.node(str(node_id), label=label)
        return dot

    ## Get the label of the node. By default, it is simply the name
    def get_dot_label(self) -> str:
        ## The name could be a full path
        name = self.cmd_invocation_with_io_vars.cmd_name
        basename = os.path.basename(str(name))
        return basename


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
        assert(False)
        if(isinstance(inputs, list)):
            self.inputs = ([], inputs)
        elif(isinstance(inputs, tuple)):
            self.inputs = inputs
        else:
            raise NotImplementedError()

    def get_input_list(self):
        inputs = self.cmd_invocation_with_io_vars.generate_inputs()
        return inputs.get_all_inputs()

    def get_output_list(self):
        return self.cmd_invocation_with_io_vars.generate_outputs()

    def get_streaming_inputs(self):
        inputs = self.cmd_invocation_with_io_vars.generate_inputs()
        return inputs.get_streaming_inputs()

    def get_configuration_inputs(self):
        inputs = self.cmd_invocation_with_io_vars.generate_inputs()
        return inputs.get_config_inputs()

    # def is_at_most_pure(self):
    #     return (self.com_category in ["stateless", "pure", "parallelizable_pure"])

    # def is_parallelizable(self):
    #     return (self.is_pure_parallelizable() or self.is_stateless())

    # def is_stateless(self):
    #     return (self.com_category == "stateless")

    # def is_pure_parallelizable(self):
    #     return (self.com_category == "parallelizable_pure")

    def is_commutative(self):
        val = self.cmd_related_properties.get_property_value('is_commutative')
        if val is not None:
            return val
        else:
            return False

    ## kk: 2021-07-23 Not totally sure if that is generally correct. Tests will say ¯\_(ツ)_/¯
    ##     I think it assumes that new options can be added in the beginning if there are no options already
    def append_options(self, new_options):
        assert(False) # unreachable
        if(len(self.com_options) > 0):
            max_opt_index = max([i for i, _opt in self.com_options])
        else:
            max_opt_index = -1
        new_com_options = [(max_opt_index + 1 + i, Arg(string_to_argument(opt))) 
                           for i, opt in enumerate(new_options)]
        self.com_options = self.com_options + new_com_options

    ## This method handles special DFG nodes specially when it has to do
    ## with turning them to commands.
    ##
    ## The goal would be for this function to be developed for more and more nodes
    ## so as to guide the second version of the annotations.
    ##
    ## TODO: Abstract this function away to annotations 2.0
    def special_to_ast(self, edges):
        assert(False) # unreachable
        # BEGIN ANNO
        return None
        # END ANNO
        ## Every argument should be completely expanded so making it a string should be fine
        if str(self.com_name) == "cat":
            redirs = self._to_ast_aux_get_redirs()
            assignments = self.com_assignments
            com_name_ast = self.com_name.to_ast()
            option_asts = [opt.to_ast() for _, opt in self.com_options]

            ## We simply turn inputs to arguments by appending them to the options
            input_arguments = self._to_ast_aux_inputs_as_args(edges, stdin_dash=True)

            ## TODO: Make sure a library of useful constructs that create
            ##       a command from a DFG node.

            ## We simply send output to stdout (as redir if needed)
            output_redir = self._to_ast_aux_single_stdout_fid(edges)

            all_arguments = [com_name_ast] + option_asts + input_arguments
            all_redirs = redirs + output_redir

            node = make_command(all_arguments, redirections=all_redirs, assignments=assignments)
            return node
        else:
            return None

    ## This function handles the input fids as arguments.
    def _to_ast_aux_inputs_as_args(self, edges, stdin_dash=False):
        assert(False) # unreachable
        input_fids = [edges[in_id][0] for in_id in self.get_input_list()]

        input_arguments = [fid.to_ast(stdin_dash=stdin_dash)
                            for fid in input_fids]
        return input_arguments

    ## This function handles the redirections when a command has a single output
    ##   and it can always be stdout.
    def _to_ast_aux_single_stdout_fid(self, edges):
        assert(False) # unreachable
        output_fids = [edges[out_id][0] for out_id in self.outputs]
        assert len(output_fids) == 1
        output_fid = output_fids[0]
        # log("output fid:", output_fid)

        output_redir = redirect_to_stdout_if_not_already(output_fid)
        # log("Redir:", output_redir)
        return output_redir

    ## Auxiliary method that returns any necessary redirections,
    ##   at the moment it doesn't look necessary.
    def _to_ast_aux_get_redirs(self):
        ## still used in to_ast
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
        return []


    ## TODO: Improve this function to be separately implemented for different special nodes,
    ##       such as cat, eager, split, etc...
    ## I do not think this separation is reasonable anymore since we remodelled nodes in a way that the back-translation is trivial
    ## One exception:
    ##  - r_wrap; currently, the wrapped command is translated at creation of the r_wrap already and
    ##    hence assumes that non-streaming inputs/outputs will not change; with a special to_ast, we could circumvent this
    def to_ast(self, edges, drain_streams):
        ## TODO: We might not want to implement this at all actually
        if (drain_streams):
            raise NotImplementedError()
        else:
            # commented since "see above"
            ## Handle special node to ast here
            # node = self.special_to_ast(edges)
            # if node is not None:
            #     return node

            redirs = self._to_ast_aux_get_redirs()
            assignments = self.com_assignments

            node = to_node_cmd_inv_with_io_vars(self.cmd_invocation_with_io_vars, edges, redirs, assignments)
            # TODO: think about redirections
            # old code for this:
            # rest_argument_fids, new_redirs = create_command_arguments_redirs(com_name_ast,
            #                                                                  option_asts,
            #                                                                  input_fids,
            #                                                                  output_fids)
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
                for i in range(len(self.get_output_list())):
                    output_edge_id = self.get_output_list()[i]
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


    ## This renames the from_id (wherever it exists in inputs or outputs)
    ## to the to_id.
    ##
    ## TODO: Make sure we don't need to change redirections here.
    ##
    ## TODO: Make this a method of graph to change the from, to too.
    def replace_edge(self, from_id, to_id):
        self.cmd_invocation_with_io_vars.replace_var(from_id, to_id)

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

    def set_used_parallelizer(self, parallelizer):
        assert(False)
        # TODO: instantiate in __init__ already in some way
        self.used_parallelizer = parallelizer

    def get_used_parallelizer(self):
        assert(False)
        return self.used_parallelizer

    def get_option_implemented_round_robin_parallelizer(self):
        for parallelizer in self.parallelizer_list:
            splitter = parallelizer.get_splitter()
            mapper_spec = parallelizer.get_mapper_spec()
            aggregator_spec = parallelizer.get_aggregator_spec()
            if splitter.is_splitter_round_robin() and mapper_spec.is_implemented and aggregator_spec.is_implemented:
                return parallelizer
        return None

    @staticmethod
    def make_simple_dfg_node_from_cmd_inv_with_io_vars(cmd_inv_with_io_vars):
        return DFGNode(cmd_inv_with_io_vars)