import copy
import annotations
from command_categories import *
from util import *
from ir_utils import *
from definitions.ir.redirection import *
from definitions.ir.resource import *

# BEGIN ANNO
from util_new_annotations import construct_property_container_from_list_of_properties
from util_new_cmd_invocations import to_node_cmd_inv_with_io_vars
from util_new_parsing import get_ast_for_flagoption, get_ast_for_argstringtype, fix_parsing_newline
from datatypes_new.BasicDatatypes import Flag, ArgStringType
# from util_new_cmd_invocations import get_command_invocation_prefix_from_dfg_node


import sys
sys.path.insert(1, "/home/felix/git-repos/MIT/annotations")
from util_new import return_empty_list_if_none_else_itself, return_default_if_none_else_itself
# END ANNO

# BEGIN REMODEL
from definitions.remodel.IOVar import IOVar
from definitions.remodel.OptionDFG import OptionDFG
from typing import List, Union, Optional

# END REMODEL

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
    ## com_options : list of tuples with the option index and the argument Arg
    ## com_redirs : list of redirections
    ## com_assignments : list of assignments
    def __init__(self,
                 cmd_invocation_with_io_vars,
                 # inputs,
                 # outputs,
                 # com_name,
                 # # com_category = None,
                 # com_options = [],
                 com_redirs = [],
                 com_assignments=[],
                 # # BEGIN ANNO
                 # flag_option_list=None,
                 # positional_config_list=None,
                 # positional_input_list=None,
                 # positional_output_list=None,
                 # implicit_use_of_stdin=False,
                 # implicit_use_of_stdout=False,
                 parallelizer_list=None,
                 cmd_related_properties=None,
                 # # END ANNO
                 # # BEGIN REMODEL
                 # stdin_used_rem : Optional[IOVar] = None,
                 # stdout_used_rem : Optional[IOVar] = None,
                 # flag_option_list_rem : List[Union[Flag, OptionDFG]] = None,
                 # operand_list_rem : List[Union[ArgStringType, IOVar]] = None
                 # # END REMODEL
                 ):
        ## Add a unique identifier to each DFGNode since id() is not guaranteed to be unique for objects that have different lifetimes.
        ## This leads to issues when nodes are deleted and new ones are created, leading to id() clashes between them
        self.id = DFGNode.next_id
        DFGNode.next_id += 1

        # self.set_inputs(inputs)
        # self.outputs = outputs
        # self.com_name = com_name
        # self.com_category = com_category
        # self.com_options = com_options # used for Mp/Agg and in to_ast
        self.com_redirs = [Redirection(redirection) for redirection in com_redirs]
        self.com_assignments = com_assignments
        # BEGIN ANNO
        # Assumption: config_list and option arguments only contains strings, i.e. of type ArgStringType
        # self.flag_option_list = return_empty_list_if_none_else_itself(flag_option_list)
        # self.positional_config_list = return_empty_list_if_none_else_itself(positional_config_list)
        # self.positional_input_list = return_empty_list_if_none_else_itself(positional_input_list)
        # self.positional_output_list = return_empty_list_if_none_else_itself(positional_output_list)
        # self.implicit_use_of_stdin = implicit_use_of_stdin
        # self.implicit_use_of_stdout = implicit_use_of_stdout
        self.parallelizer_list = return_empty_list_if_none_else_itself(parallelizer_list)
        default_cmd_properties = construct_property_container_from_list_of_properties([])
        self.cmd_related_properties = return_default_if_none_else_itself(cmd_related_properties, default_cmd_properties)
        self.cmd_invocation_with_io_vars = cmd_invocation_with_io_vars
        # END ANNO

        # log("Node created:", self.id, self)

    def __repr__(self):
        ## BEGIN ANNO
        # NEW
        # TODO ANNO
        return str(self.cmd_invocation_with_io_vars)
        # OLD
        # # prefix = "Node"
        # # if (self.com_category == "stateless"):
        # #     prefix = "Stateless"
        # # elif (self.com_category == "pure"):
        # #     prefix = "Pure"
        # # elif (self.is_pure_parallelizable()):
        # #     prefix = "Par. Pure"
        # # if (self.is_commutative()):
        # #     prefix = 'Commutative ' + prefix
        # # output = "{}: \"{}\" in:{} out:{}".format(
        # #     prefix, self.com_name,
        # #     self.get_input_list(),
        # #     self.outputs)
        # return output
        ## END ANNO

    ## Generates a dot node for the DFG node
    def add_dot_node(self, dot, node_id):
        label = self.get_dot_label()
        dot.node(str(node_id), label=label)
        return dot

    ## Get the label of the node. By default, it is simply the name
    def get_dot_label(self) -> str:
        ## The name could be a full path
        name = self.com_name
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

    ## BEGIN ANNO
    # def is_at_most_pure(self):
    #     return (self.com_category in ["stateless", "pure", "parallelizable_pure"])

    def is_parallelizable(self):
        return (self.is_pure_parallelizable() or self.is_stateless())

    def is_stateless(self):
        return (self.com_category == "stateless")

    def is_pure_parallelizable(self):
        return (self.com_category == "parallelizable_pure")
    # END ANNO

    def is_commutative(self):
        # BEGIN ANNO
        # OLD
        # return ('commutative' in self.com_properties)
        # NEW
        val = self.cmd_related_properties.get_property_value('commutative')
        if val is not None:
            return val
        else:
            return False
        # END ANNO

    ## kk: 2021-07-23 Not totally sure if that is generally correct. Tests will say ¯\_(ツ)_/¯
    ##     I think it assumes that new options can be added in the beginning if there are no options already
    def append_options(self, new_options):
        assert(False)
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
        input_fids = [edges[in_id][0] for in_id in self.get_input_list()]

        input_arguments = [fid.to_ast(stdin_dash=stdin_dash)
                            for fid in input_fids]
        return input_arguments

    ## This function handles the redirections when a command has a single output
    ##   and it can always be stdout.
    def _to_ast_aux_single_stdout_fid(self, edges):
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


    ## TODO: Improve this functio to be separately implemented for different special nodes,
    ##       such as cat, eager, split, etc...
    def to_ast(self, edges, drain_streams):    
        ## TODO: We might not want to implement this at all actually
        if (drain_streams):
            raise NotImplementedError()
        else:
            ## Handle special node to ast here
            node = self.special_to_ast(edges)
            if node is not None:
                return node
            

            redirs = self._to_ast_aux_get_redirs()
            assignments = self.com_assignments
            ## Start filling in the arguments
            # opt_arguments = []
            # BEGIN ANNO
            # LOG
            # log(f'com_name: {self.com_name}')
            # log(f'edges: {edges}')
            # log(f'inputs: {self.inputs}')
            # log(f'outputs: {self.outputs}')
            # log(f'com_redirs: {self.com_redirs}')
            # log(f'pos config: {self.positional_config_list}')
            # log(f'pos input: {self.positional_input_list}')
            # log(f'pos output: {self.positional_output_list}')
            # log(f'com_options: {self.com_options}')
            # log(f'flag_option_list: {self.flag_option_list}')
            #OLD
            # for i, opt in self.com_options:
            #     ## Pad the argument list with None
            #     opt_arguments = pad(opt_arguments, i)
            #     opt_arguments[i] = opt.to_ast()
            # log(f'opt_arguments: {format_args([val for val in opt_arguments if val is not None])}')
            # NEW
            # TODO: FIX this?!
            # opt_arguments = [get_ast_for_flagoption(flagoption) for flagoption in self.flag_option_list]
            # positional_config_fixed_parsing_newline = [fix_parsing_newline(arg) for arg in self.positional_config_list]
            # opt_arguments += [get_ast_for_argstringtype(arg) for arg in positional_config_fixed_parsing_newline]
            # log(f'opt_arguments_new: {format_args(opt_arguments)}')
            # END ANNO

            node = to_node_cmd_inv_with_io_vars(self.cmd_invocation_with_io_vars, edges, redirs, assignments)
            # TODO: think about redirections

            # com_name_ast = self.com_name.to_ast()

            ##
            ## 1. Find the input and output fids
            ## 2. Construct the rest of the arguments and input/output redirections according to
            ##    the command IO
            # input_fids = [edges[in_id][0] for in_id in self.get_input_list()]
            # output_fids = [edges[out_id][0] for out_id in self.outputs]
            ## FS: com_option is still used but with new model for nodes, this function shall go
            # option_asts = [opt.to_ast() for _, opt in self.com_options]
            # rest_argument_fids, new_redirs = create_command_arguments_redirs(com_name_ast,
            #                                                                  option_asts,
            #                                                                  input_fids,
            #                                                                  output_fids)
            
            ## Transform the rest of the argument fids to arguments
            ## Since some of the rest_arguments can be None (they only contain inputs and outputs)
            ## we need to make sure that we don't turn None objects to asts.
            ##
            ## The None fields need to be filtered out because they are taken care of by the interleave function.
            ##
            ## TODO: Is this actually OK?
            # rest_arguments = [fid.to_ast()
            #                   for fid in rest_argument_fids
            #                   if not fid is None]
            # log(f'rest_arguments: {format_args(rest_arguments)}')

            ## Interleave the arguments since options args might contain gaps.
            # BEGIN ANNO
            # OLD
            # arguments = interleave_args(opt_arguments, rest_arguments)
            # NEW
            # arguments = opt_arguments + rest_arguments
            # log(f'arguments: {format_args(arguments)}')
            # END ANNO

            # all_arguments = [com_name_ast] + arguments
            # log(f'all arguments: {format_args(all_arguments)}')
            # log(f'\n\n')
            # all_redirs = redirs + new_redirs

            # node = make_command(all_arguments, redirections=all_redirs, assignments=assignments)
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


    ## This renames the from_id (wherever it exists in inputs or outputs)
    ## to the to_id.
    ##
    ## TODO: Make sure we don't need to change redirections here.
    ##
    ## TODO: Make this a method of graph to change the from, to too.
    def replace_edge(self, from_id, to_id):
        self.cmd_invocation_with_io_vars.replace_var(from_id, to_id)
        # new_config_inputs = self.replace_edge_in_list(self.inputs[0], from_id, to_id)
        # new_standard_inputs = self.replace_edge_in_list(self.inputs[1], from_id, to_id)
        # new_outputs = self.replace_edge_in_list(self.outputs, from_id, to_id)

        # self.set_inputs((new_config_inputs, new_standard_inputs))
        # self.outputs = new_outputs

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
        # TODO: instantiate in __init__ already in some way
        self.used_parallelizer = parallelizer

    def get_used_parallelizer(self):
        return self.used_parallelizer

    def get_option_round_robin_parallelizer(self):
        for parallelizer in self.parallelizer_list:
            splitter = parallelizer.get_splitter()
            if splitter.is_splitter_round_robin():
                return parallelizer
        return None

    @staticmethod
    def make_simple_dfg_node_from_cmd_inv_with_io_vars(cmd_inv_with_io_vars):
        return DFGNode(cmd_inv_with_io_vars)