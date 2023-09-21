import sys

from config import get_path_annotation_repo
sys.path.insert(1, get_path_annotation_repo())
from datatypes_new.CommandInvocationInitial import CommandInvocationInitial
from datatypes_new.BasicDatatypes import ArgStringType
from datatypes_new.BasicDatatypesWithIO import FileNameWithIOInfo, StdDescriptorWithIOInfo, OptionWithIO
from annotation_generation_new.datatypes.InputOutputInfo import InputOutputInfo
from annotation_generation_new.datatypes.ParallelizabilityInfo import ParallelizabilityInfo
from datatypes_new.CommandInvocationWithIOVars import CommandInvocationWithIOVars

from annotations_utils.util_parsing import parse_arg_list_to_command_invocation
from annotations_utils.util_cmd_invocations import get_input_output_info_from_cmd_invocation_util, get_parallelizability_info_from_cmd_invocation_util
from annotations_utils.util_mapper import get_mapper_as_dfg_node_from_node, get_map_output_files
from annotations_utils.util_aggregator import get_aggregator_as_dfg_node_from_node
from annotations_utils.util_file_descriptors import resource_from_file_descriptor

from definitions.ir.file_id import *
from definitions.ir.nodes.cat import *

import definitions.ir.nodes.pash_split as pash_split
import definitions.ir.nodes.r_merge as r_merge
import definitions.ir.nodes.r_split as r_split
import definitions.ir.nodes.r_wrap as r_wrap
import definitions.ir.nodes.r_unwrap as r_unwrap

from ir_utils import *
from util import *

import config

## Creates a file id for a given resource
def create_file_id_for_resource(resource, fileIdGen):
    file_id = create_split_file_id(fileIdGen)
    file_id.set_resource(resource)
    return file_id

## Creates a file id that has a given maximum length
def create_split_file_id(fileIdGen):
    file_id = fileIdGen.next_file_id()
    return file_id

class FileIdGen:
    def __init__(self, next = 0, prefix = ""):
        self.next = next + 1
        directory = f"{str(uuid.uuid4().hex)}"
        self.prefix = f"{directory}/{prefix}"
        directory_path = os.path.join(config.PASH_TMP_PREFIX, self.prefix)  
        os.makedirs(directory_path)

    def next_file_id(self):
        fileId = FileId(self.next, self.prefix)
        self.next += 1
        return fileId

    def next_temporary_file_id(self):
        fileId = self.next_file_id()
        fileId.make_temporary_file()
        return fileId

    def next_ephemeral_file_id(self):
        fileId = self.next_file_id()
        fileId.make_ephemeral()
        return fileId

    def bump_counter_to_value_of(self, OtherFileIdGen):
        # TODO: find a better solution to make unique numbers, currently: set to max-value + 1
        self.next = OtherFileIdGen.next + 1

## Returns the resource or file descriptor related to this specific opt_or_fd
## NOTE: Assumes that everything is expanded. 
def get_option_or_fd(opt_or_fd, options, fileIdGen):
    if(isinstance(opt_or_fd, tuple)
       and len(opt_or_fd) == 2
       and opt_or_fd[0] == "option"):
        resource = FileResource(Arg(options[opt_or_fd[1]]))
    else:
        ## TODO: Make this be a subtype of Resource
        if(opt_or_fd == "stdin"):
            resource = ("fd", 0)
        elif(opt_or_fd == "stdout"):
            resource = ("fd", 1)
        elif(opt_or_fd == "stderr"):
            resource = ("fd", 2)
        else:
            raise NotImplementedError()
        resource = FileDescriptorResource(resource)
    
    fid = create_file_id_for_resource(resource, fileIdGen)
    return fid

## Get the options as arguments
def get_option(opt_or_fd, options, fileIdGen):
    assert(isinstance(opt_or_fd, tuple)
       and len(opt_or_fd) == 2
       and opt_or_fd[0] == "option")
    arg = Arg(options[opt_or_fd[1]])
    return (opt_or_fd[1], arg)

## This function 
def create_edges_from_opt_or_fd_list(opt_or_fd_list, edges_dict, options, fileIdGen):
    new_edge_list = []
    for opt_or_fd in opt_or_fd_list:
        fid = get_option_or_fd(opt_or_fd, options, fileIdGen)
        fid_id = fid.get_ident()
        edges_dict[fid_id] = (fid, None, None)
        new_edge_list.append(fid_id)
    return new_edge_list


def find_input_edges(positional_input_list, implicit_use_of_stdin, dfg_edges, fileIdGen) -> List[int]:
    assert (not implicit_use_of_stdin or len(positional_input_list) == 0)
    if implicit_use_of_stdin:
        resources = [FileDescriptorResource(("fd", 0))]
    else:
        resources = [resource_from_file_descriptor(input_el) for input_el in positional_input_list]
    file_ids = [create_file_id_for_resource(resource, fileIdGen) for resource in resources]
    return get_edge_list_from_file_id_list(dfg_edges, file_ids)


def find_output_edges(positional_output_list, implicit_use_of_stdout, dfg_edges, fileIdGen) -> List[int]:
    assert (not implicit_use_of_stdout or len(positional_output_list) == 0)
    if implicit_use_of_stdout:
        resources = [FileDescriptorResource(("fd", 1))]
    else:
        resources = [resource_from_file_descriptor(input_el) for input_el in positional_output_list]
    file_ids = [create_file_id_for_resource(resource, fileIdGen) for resource in resources]
    return get_edge_list_from_file_id_list(dfg_edges, file_ids)


def get_edge_list_from_file_id_list(dfg_edges, file_ids):
    new_edge_list = []
    for file_id in file_ids:
        fid_id = file_id.get_ident()
        dfg_edges[fid_id] = (file_id, None, None)
        new_edge_list.append(fid_id)
    return new_edge_list


def add_file_id_vars(command_invocation_with_io, fileIdGen):
    # make pass over everything and create file_id for everything
    # only for operands for now:
    dfg_edges = {}
    new_flagoption_list = []
    new_operand_list = []
    access_map = dict()

    def add_var_for_descriptor(operand):
        resource = resource_from_file_descriptor(operand)
        file_id = create_file_id_for_resource(resource, fileIdGen)
        fid_id = file_id.get_ident()
        dfg_edges[fid_id] = (file_id, None, None)
        access_map[fid_id] = operand.get_access()
        return fid_id

    for i in range(len(command_invocation_with_io.flag_option_list)):
        flagoption = command_invocation_with_io.flag_option_list[i]
        if isinstance(flagoption, OptionWithIO) and not isinstance(flagoption.option_arg, ArgStringType):
            fid_id = add_var_for_descriptor(flagoption.option_arg)
            new_option = OptionWithIOVar(flagoption.name, fid_id)
            new_flagoption_list.append(new_option)
        else: # Flag
            new_flagoption_list.append(flagoption)

    for i in range(len(command_invocation_with_io.operand_list)):
        operand = command_invocation_with_io.operand_list[i]
        if isinstance(operand, FileNameWithIOInfo) or isinstance(operand, StdDescriptorWithIOInfo):
            fid_id = add_var_for_descriptor(operand)
            new_operand_list.append(fid_id)
        else:
            new_operand_list.append(operand)
    if command_invocation_with_io.implicit_use_of_streaming_input:
        new_implicit_use_of_streaming_input = add_var_for_descriptor(command_invocation_with_io.implicit_use_of_streaming_input)
    else:
        new_implicit_use_of_streaming_input = None
    if command_invocation_with_io.implicit_use_of_streaming_output:
        new_implicit_use_of_streaming_output = add_var_for_descriptor(command_invocation_with_io.implicit_use_of_streaming_output)
    else:
        new_implicit_use_of_streaming_output = None

    command_invocation_with_io_vars = CommandInvocationWithIOVars(cmd_name=command_invocation_with_io.cmd_name,
                                       flag_option_list=new_flagoption_list,
                                       operand_list=new_operand_list,
                                       implicit_use_of_streaming_input=new_implicit_use_of_streaming_input,
                                       implicit_use_of_streaming_output=new_implicit_use_of_streaming_output,
                                       access_map=access_map)
    return command_invocation_with_io_vars, dfg_edges


def compile_command_to_DFG(fileIdGen, command, options,
                           redirections=[]):
    command_invocation: CommandInvocationInitial = parse_arg_list_to_command_invocation(command, options)
    io_info: InputOutputInfo = get_input_output_info_from_cmd_invocation_util(command_invocation)
    para_info: ParallelizabilityInfo = get_parallelizability_info_from_cmd_invocation_util(command_invocation)
    if para_info is None:
        para_info = ParallelizabilityInfo() # defaults to no parallelizer's and all properties False
    command_invocation_with_io = io_info.apply_input_output_info_to_command_invocation(command_invocation)
    if para_info is None:
        para_info = ParallelizabilityInfo()  # defaults to no parallelizer's and all properties False
    parallelizer_list, round_robin_compatible_with_cat, is_commutative = para_info.unpack_info()
    property_list = [('round_robin_compatible_with_cat', round_robin_compatible_with_cat),
                     ('is_commutative', is_commutative)]
    cmd_related_properties = construct_property_container_from_list_of_properties(property_list)

    ## TODO: Make an empty IR and add edges and nodes incrementally (using the methods defined in IR).

    ## Add all inputs and outputs to the DFG edges
    cmd_invocation_with_io_vars, dfg_edges = add_file_id_vars(command_invocation_with_io, fileIdGen)
    com_redirs = redirections
    ## TODO: Add assignments
    com_assignments = []

    ## TODO: Combine them both in a constructor that decided whether to instantiate Cat or DFGNode
    # if(str(com_name) == "cat"):
    #     dfg_node = Cat(dfg_inputs,
    #                    dfg_outputs,
    #                    com_name,
    #                    ## TODO: We don't really need to pass category, name, or input_consumption for Cat
    #                    com_category,
    #                    com_options=dfg_options,
    #                    com_redirs=com_redirs,
    #                    com_assignments=com_assignments,
    #                    )
    # elif(str(com_name) == "hdfs" and str(dfg_options[0][1]) == "dfs" and str(dfg_options[1][1]) == "-cat"):
    #     dfg_node = HDFSCat(dfg_inputs,
    #                     dfg_outputs,
    #                     com_name,
    #                     com_category,
    #                     com_options=dfg_options,
    #                     com_redirs=com_redirs,
    #                     com_assignments=com_assignments)
    # else:
    if(True):
        ## Assume: Everything must be completely expanded
        ## TODO: Add an assertion about that.
        dfg_node = DFGNode(cmd_invocation_with_io_vars,
                           com_redirs=com_redirs,
                           com_assignments=com_assignments,
                           parallelizer_list=parallelizer_list,
                           cmd_related_properties=cmd_related_properties
                           )

    # if(not dfg_node.is_at_most_pure()): # which consequences has this check had?
    #     raise ValueError()

    node_id = dfg_node.get_id()

    ## Assign the from, to node in edges
    for fid_id in dfg_node.get_input_list():
        fid, from_node, to_node = dfg_edges[fid_id]
        assert(to_node is None)
        dfg_edges[fid_id] = (fid, from_node, node_id)
    
    for fid_id in dfg_node.get_output_list():
        fid, from_node, to_node = dfg_edges[fid_id]
        assert(from_node is None)
        dfg_edges[fid_id] = (fid, node_id, to_node)
    
    dfg_nodes = {node_id : dfg_node}
    dfg = IR(dfg_nodes, dfg_edges)
    return dfg


def make_split_files(input_id, fan_out, fileIdGen, r_split_flag, r_split_batch_size):
    assert(fan_out > 1)
    ## Generate the split file ids
    out_fids = [fileIdGen.next_file_id() for i in range(fan_out)]
    out_ids = [fid.get_ident() for fid in out_fids]
    split_com = make_split_file(input_id, out_ids, r_split_flag, r_split_batch_size)
    return [split_com], out_fids

def make_split_file(input_id, out_ids, r_split_flag, r_split_batch_size):
    if(r_split_flag):
        split_com = r_split.make_r_split(input_id, out_ids, r_split_batch_size)
    else:
        split_com = pash_split.make_split_file(input_id, out_ids)
    return split_com

##
## Node builder functions
##

def make_tee(input, outputs):
    com_name = Arg(string_to_argument("tee"))
    com_category = "pure"
    return DFGNode([input],
                   outputs,
                   com_name, 
                   com_category)

def make_map_node(node, new_inputs, new_outputs, parallelizer):
    return get_mapper_as_dfg_node_from_node(node, parallelizer, new_inputs, new_outputs)

## Makes a wrap node that encloses a map parallel node.
##
## At the moment it only works with one input and one output since wrap cannot redirect input in the command.
def make_wrap_map_node(node, new_inputs, new_outputs):
    assert(len(new_inputs) == 1)
    assert(len(new_outputs) == 1)

    new_node = make_map_node(node, new_inputs, new_outputs)
    wrap_node = r_wrap.wrap_node(new_node)
    return wrap_node



## Note: This might need more information. E.g. all the file
## descriptors of the IR, and in general any other local information
## that might be relevant.
class IR:

    ## TODO: Embed the fileIdGen as a field of the IR

    ## IR Assumptions:
    ##
    ## - Each node has a list of incoming files in order of
    ##   consumption.
    ##
    ## - If two nodes have the same file as output, then they both
    ##   write to it concurrently.
    def __init__(self, nodes, edges, background = False):
        self.nodes = nodes
        self.edges = edges
        self.background = background

        ## Apply the redirections for each separate node.
        ## This needs to be called here because nodes do not
        ## have information about the edges on their own.
        self.apply_redirections()

    def __repr__(self):
        output = "(|-{} IR: {} {}-|)".format(self.get_stdin(), list(self.nodes.values()), self.get_stdout())
        return output

    ## Initialize all edges
    def apply_redirections(self):
        for _, node in self.nodes.items():
            node.apply_redirections(self.edges)
        
        ## We need to merge common files after redirections have been applied.
        self.combine_common_files()

    ## Refactor these to call .add_edge, and .set_edge_to/from 
    ## Add an edge that points to a node
    def add_to_edge(self, to_edge, node_id):
        edge_id = to_edge.get_ident()
        assert(not edge_id in self.edges)
        self.edges[edge_id] = (to_edge, None, node_id)

    ## Add an edge that starts from a node
    def add_from_edge(self, node_id, from_edge):
        edge_id = from_edge.get_ident()
        assert(not edge_id in self.edges)
        self.edges[edge_id] = (from_edge, node_id, None)

    def set_edge_to(self, edge_id, to_node_id):
        edge_fid, from_node, old_to_node = self.edges[edge_id]
        self.edges[edge_id] = (edge_fid, from_node, to_node_id)

    def set_edge_from(self, edge_id, from_node_id):
        edge_fid, old_from_node, to_node = self.edges[edge_id]
        self.edges[edge_id] = (edge_fid, from_node_id, to_node)

    def get_edge_fid(self, fid_id):
        if(fid_id in self.edges):
            return self.edges[fid_id][0]
        else:
            return None

    def get_edge_from(self, edge_id):
        if(edge_id in self.edges):
            return self.edges[edge_id][1]
        else:
            return None

    def replace_edge(self, old_edge_id, new_edge_fid):
        assert(new_edge_fid not in self.all_fids())
        new_edge_id = new_edge_fid.get_ident()
        old_fid, from_node, to_node = self.edges[old_edge_id]
        self.edges[new_edge_id] = (new_edge_fid, from_node, to_node)
        if from_node:
            self.get_node(from_node).replace_edge(old_edge_id, new_edge_id)
        if to_node:
            self.get_node(to_node).replace_edge(old_edge_id, new_edge_id)
        del self.edges[old_edge_id]
        
    def get_stdin(self):
        stdin_id = self.get_stdin_id()
        stdin_fid = self.get_edge_fid(stdin_id)
        return stdin_fid

    def get_stdout(self):
        stdout_id = self.get_stdout_id()
        stdout_fid = self.get_edge_fid(stdout_id)
        return stdout_fid

    ## Gets the fid that points to the stdin of this DFG
    def get_stdin_id(self):
        ## ASSERT: There must be only one
        stdin_id = None
        for edge_id, (edge_fid, _from, _to) in self.edges.items():
            resource = edge_fid.get_resource()
            if(resource.is_stdin()):
                assert(stdin_id is None)
                stdin_id = edge_id
        return stdin_id  

    def get_stdout_id(self):
        ## ASSERT: There must be only one
        stdout_id = None
        for edge_id, (edge_fid, _from, _to) in self.edges.items():
            resource = edge_fid.get_resource()
            if(resource.is_stdout()):
                # This is not true when using distributed_exec
                # assert(stdout_id is None)
                stdout_id = edge_id
        return stdout_id   

    def serialize(self):
        output = "Nodes:\n"
        all_file_ids = ""
        for i, node in enumerate(self.nodes):
            serialized_input_file_ids = " ".join([fid.serialize()
                                                  for fid in node.get_input_file_ids()])
            serialized_output_file_ids = " ".join([fid.serialize()
                                                   for fid in node.get_output_file_ids()])
            all_file_ids += serialized_input_file_ids + " "
            all_file_ids += serialized_output_file_ids + " "
            output += "{} in: {} out: {} command: {}\n".format(i, serialized_input_file_ids,
                                                               serialized_output_file_ids,
                                                               node.serialize())
        output = "File ids:\n{}\n".format(all_file_ids) + output
        return output


    def to_ast(self, drain_streams):
        asts = []

        ## Initialize the pids_to_kill variable
        asts.append(self.init_pids_to_kill())

        fileIdGen = self.get_file_id_gen()

        ## Redirect stdin
        stdin_id = self.get_stdin_id()
        if (not stdin_id is None):
            ## Create a new ephemeral resource to redirect stdin to.
            fid = fileIdGen.next_file_id()
            fid.make_ephemeral()
            file_to_redirect_to = fid.to_ast()
            ## Change the stdin_id to point to this resource
            _prev_fid, from_node, to_node = self.edges[stdin_id]
            self.edges[stdin_id] = (fid, from_node, to_node)
            ## Create a command that redirects stdin to this ephemeral fid
            redirect_stdin_script = os.path.join(config.PASH_TOP, config.config['runtime']['redirect_stdin_binary'])
            com_args = [string_to_argument('source'), string_to_argument(redirect_stdin_script), file_to_redirect_to]
            com = make_command(com_args)
            asts.append(com)

        ## Make the dataflow graph
        ##
        ## TODO: Normally this should have all sink nodes at the end, but
        ##       for now we just have the stdout node in the end 
        ##       (since this is always the output in our benchmarks).
        # sink_node_ids = self.sink_nodes()
        ##
        ## TODO: Support more than one output (and less than one).
        ##       For this we need to update wait.
        ##
        ## For now we just allow more than one output by waiting for one of them
        ## at random.
        stdout_edge_id = self.get_stdout_id()
        if (not stdout_edge_id is None):
            sink_node_ids = [self.edges[stdout_edge_id][1]]
        else:
            sink_node_ids = self.sink_nodes()
            sink_node_ids = [sink_node_ids[0]]


        for node_id, node in self.nodes.items():
            if(not node_id in sink_node_ids):
                node_ast = node.to_ast(self.edges, drain_streams)
                asts.append(make_background(node_ast))
                ## Gather all pids
                assignment = self.collect_pid_assignment()
                asts.append(assignment)

        ## Put the output node in the end for wait to work.
        for node_id in sink_node_ids:
            node = self.get_node(node_id)
            node_ast = node.to_ast(self.edges, drain_streams)
            asts.append(make_background(node_ast))
            ## Gather all pids
            assignment = self.collect_pid_assignment()
            asts.append(assignment)

        return asts
    
    def collect_pid_assignment(self):
        ## Creates:
        ## pids_to_kill="$! $pids_to_kill"
        var_name = 'pids_to_kill'
        rval = quote_arg([standard_var_ast('!'),
                          char_to_arg_char(' '),
                          standard_var_ast(var_name)])
        return make_assignment(var_name, [rval])
    
    def init_pids_to_kill(self):
        ## Creates:
        ## pids_to_kill=""
        var_name = 'pids_to_kill'
        rval = quote_arg([])
        return make_assignment(var_name, [rval])

    ## TODO: Delete this
    def set_ast(self, ast):
        self.ast = ast

    def set_background(self, background):
        self.background = background

        if (background):
            ## Since the IR is in the background, we don't have access to
            ## its stdin, stdout anymore
            self.stdin = []
            self.stdout = []

    def is_in_background(self):
        return self.background

    def pipe_append(self, other):
        assert(self.valid())
        assert(other.valid())

        ## This combines the two IRs by adding all of the nodes
        ## together, and by union-ing the stdout of the first with the
        ## stdin of the second.
        ##
        ## Question: What happens if one of them is NULL. This
        ##           shouldn't be the case after we check that
        ##           both self and other are not empty.
        my_out = self.get_stdout_id()
        other_in = other.get_stdin_id()
        assert(not my_out is None)
        assert(not other_in is None)


        _other_in_fid, from_node, other_in_node_id = other.edges[other_in]
        assert(from_node is None)
        ## ... = OtherInNode(..., other_in, ...)
        ##          v
        ## ... = OtherInNode(..., my_out, ...)
        other_in_node = other.nodes[other_in_node_id]
        other_in_node.replace_edge(other_in, my_out)
        other.edges.pop(other_in)

        ## Make the my_out id to be ephemeral file.
        my_out_fid, from_node, to_node = self.edges[my_out]
        assert(to_node is None)
        my_out_fid.make_ephemeral()

        ## Add the other node in my edges
        self.edges[my_out] = (my_out_fid, from_node, other_in_node_id)

        ## Just call union here
        self.union(other)

    def background_union(self, other):
        assert(self.valid())
        assert(other.valid())
        assert(self.is_in_background())
        ## This combines two IRs where at least the first one is in
        ## background. This means that the stdin only works with the second
        ## the second (or None if both are in background). Also if
        ## both are in background, their union is also in background.

        ## If one of them is not in the background, then the whole
        ## thing isn't.
        if (not other.is_in_background()):
            self.set_background(other.is_in_background())

        self.union(other)

    def union(self, other):
        ## Merge the nodes of the two DFGs
        all_nodes = {**self.nodes, **other.nodes}

        ## Merge edges
        all_edges = {**self.edges, **other.edges}

        ## TODO: Check that all ids are OK (no cycles etc)
        self.nodes = all_nodes
        self.edges = all_edges

        ## TODO: Handle connections of common files (pipes, etc)
        self.combine_common_files()


    ## Combines (unions) files that refer to the same resource.
    ##
    ## WARNING: This assumes that comparing file names statically
    ## (syntactically) for equality, implies semantic
    ## equality. However this is not true. There are cases where
    ## different identifiers could refer to the same file.
    ##
    ## TODO: Fix the above issue by ensuring normalized absolute names
    ##
    ## Q: Are there also cases where a same name (let's say a
    ## variable) could point to different files in different parts of
    ## the IR? Maybe it can be true if a command is run with
    ## variable assignments)
    def combine_common_files(self):

        ## For now we just unify a file if it exists exactly twice,
        ## once at the input of a node and once at the output of
        ## another node. If a file exists in several input locations,
        ## we don't unify it. Also if a file exists in more than 1
        ## input and 1 output (or more than 1 output in general) we
        ## signal an error.

        ## For all inputs of all nodes, check if they are the output
        ## of exactly one other node.
        # log("Combining files for:", self)
        for node_id1, _node1 in self.nodes.items():
            inputs_with_file_resource = [(id, fid) for id, fid in self.get_node_input_ids_fids(node_id1)
                                         if fid.has_file_resource()]
            for id_in, fid_in in inputs_with_file_resource:
                in_resource = fid_in.get_resource()
                number_of_out_resources = 0
                for node_id2, _node2 in self.nodes.items():
                    outputs_with_file_resource = [(id, fid) for id, fid in self.get_node_output_ids_fids(node_id2)
                                                  if fid.has_file_resource()]
                    for id_out, fid_out in outputs_with_file_resource:
                        out_resource = fid_out.get_resource()
                        ## Do not combine if the ids of the edges are already the same
                        if (not id_in == id_out
                            and in_resource == out_resource):
                            number_of_out_resources += 1
                            ## They point to the same File resource so we need to unify their fids
                            self.nodes[node_id2].replace_edge(id_out, id_in)
                            self.set_edge_from(id_in, node_id2)
                            self.set_edge_from(id_out, None)

                ## Exit with an error if a file is written by more than one node.
                ##
                ## TODO: Could this ever be improved for additional performance?
                assert(number_of_out_resources <= 1)

    ## Returns all the file identifiers in the IR.
    def all_fids(self):
        all_fids = [fid for fid, _from_node, _to_node in self.edges.values()]
        return all_fids

    ## Returns all input fids of the IR
    def all_input_fids(self):
        all_input_fids = [fid for fid, from_node, _to_node in self.edges.values()
                          if from_node is None]
        return all_input_fids

    ## Returns all output fids of the IR
    def all_output_fids(self):
        all_output_fids = [fid for fid, _from_node, to_node in self.edges.values()
                          if to_node is None]
        return all_output_fids

    ## Returns the sources of the IR.
    ##   This includes both the nodes that have an incoming edge (file) that has no from_node,
    ##     but also nodes that have no incoming edge (generator nodes). 
    def source_nodes(self):
        sources = set()
        for _edge_fid, from_node, to_node in self.edges.values():
            if(from_node is None and not to_node is None):
                sources.add(to_node)
        for node_id, node in self.nodes.items():
            if len(node.get_input_list()) == 0:
                sources.add(node_id)
        return list(sources)

    def sink_nodes(self):
        sources = set()
        for _edge_fid, from_node, to_node in self.edges.values():
            if(to_node is None and not from_node is None):
                sources.add(from_node)
        return list(sources)

    def get_node_inputs(self, node_id):
        input_edge_ids = self.nodes[node_id].get_input_list()
        return input_edge_ids

    def get_node_outputs(self, node_id):
        output_edge_ids = self.nodes[node_id].get_output_list()
        return output_edge_ids

    def get_next_nodes(self, node_id):
        output_edge_ids = self.get_node_outputs(node_id)
        next_nodes = []
        for edge_id in output_edge_ids:
            _fid, from_node, to_node = self.edges[edge_id]
            assert(from_node == node_id)
            if(not to_node is None):
                next_nodes.append(to_node)
        return next_nodes

    def get_previous_nodes(self, node_id):
        input_edge_ids = self.get_node_inputs(node_id)
        previous_nodes = []
        for edge_id in input_edge_ids:
            _fid, from_node, to_node = self.edges[edge_id]
            assert(to_node == node_id)
            if(not from_node is None):
                previous_nodes.append(from_node)
        return previous_nodes

    def get_node_input_ids_fids(self, node_id):
        node = self.get_node(node_id)
        return [(input_edge_id, self.edges[input_edge_id][0]) for input_edge_id in node.get_input_list()]

    def get_node_input_ids(self, node_id):
        return [fid_id for fid_id, _fid in self.get_node_input_ids_fids(node_id)]

    def get_node_input_fids(self, node_id):
        return [fid for _fid_id, fid in self.get_node_input_ids_fids(node_id)]

    def get_node_output_ids_fids(self, node_id):
        node = self.get_node(node_id)
        return [(output_edge_id, self.edges[output_edge_id][0]) for output_edge_id in node.get_output_list()]

    def get_node_output_ids(self, node_id):
        return [fid_id for fid_id, _fid in self.get_node_output_ids_fids(node_id)]

    def get_node_output_fids(self, node_id):
        return [fid for _fid_id, fid in self.get_node_output_ids_fids(node_id)]

    def get_node(self, node_id):
        return self.nodes[node_id]

    ## This command gets all file identifiers of the graph, and
    ## returns a fileId generator that won't clash with the existing
    ## ones.
    def get_file_id_gen(self):
        max_id = max(self.edges.keys())
        return FileIdGen(max_id)

    def remove_node(self, node_id):
        node = self.nodes.pop(node_id)
        ## Remove the node in the edges dictionary
        for in_id in node.get_input_list():
            self.set_edge_to(in_id, None)

        for out_id in node.get_output_list():
            self.set_edge_from(out_id, None)


    def add_node(self, node):
        node_id = node.get_id()
        self.nodes[node_id] = node
        ## Add the node in the edges dictionary
        for in_id in node.get_input_list():
            self.set_edge_to(in_id, node_id)

        for out_id in node.get_output_list():
            self.set_edge_from(out_id, node_id)

    def generate_ephemeral_edges(self, fileIdGen, num_of_edges):
        file_ids = [fileIdGen.next_ephemeral_file_id() for _ in range(num_of_edges)]
        self.add_edges(file_ids)
        return [edge_fid.get_ident() for edge_fid in file_ids]

    def add_edges(self, edge_fids):
        for edge_fid in edge_fids:
            self.add_edge(edge_fid)

    def add_edge(self, edge_fid):
        fid_id = edge_fid.get_ident()
        assert(not fid_id in self.edges)
        self.edges[fid_id] = (edge_fid, None, None)

    ## Note: We assume that the lack of nodes is an adequate condition
    ##       to check emptiness.
    def empty(self):
        return (len(self.nodes) == 0)

    def apply_parallelization_to_node(self, node_id, parallelizer, fileIdGen, fan_out,
                                            batch_size, no_cat_split_vanish, r_split_batch_size):
        splitter = parallelizer.get_splitter()
        if splitter.is_splitter_round_robin():
            # TODO: for both functions, check which parameters are needed
            self.apply_round_robin_parallelization_to_node(node_id, parallelizer, fileIdGen, fan_out,
                                            batch_size, no_cat_split_vanish, r_split_batch_size)
        elif splitter.is_splitter_round_robin_with_unwrap_flag():
            self.apply_round_robin_with_unwrap_flag_parallelization_to_node(node_id, parallelizer, fileIdGen, fan_out,
                                                           batch_size, no_cat_split_vanish, r_split_batch_size)
        elif splitter.is_splitter_consec_chunks():
            self.apply_consecutive_chunks_parallelization_to_node(node_id, parallelizer, fileIdGen, fan_out,
                                                                 batch_size, no_cat_split_vanish, r_split_batch_size)
        else:
            raise Exception("Splitter not yet implemented")

    def apply_round_robin_parallelization_to_node(self, node_id, parallelizer, fileIdGen, fan_out,
                                                        batch_size, no_cat_split_vanish, r_split_batch_size):
        # TODO: this control flow should move done to aggregators once we implement them;
        #  currently, this cannot be done since splitter etc. would be added...
        aggregator_spec = parallelizer.get_aggregator_spec()
        if aggregator_spec.is_aggregator_spec_adj_lines_merge():
            raise Exception("adj_lines_merge not yet implemented in PaSh")
        elif aggregator_spec.is_aggregator_spec_adj_lines_seq():
            raise Exception("adj_lines_seq not yet implemented in PaSh")
        elif aggregator_spec.is_aggregator_spec_adj_lines_func():
            raise Exception("adj_lines_func not yet implemented in PaSh")
        # END of what to move

        node = self.get_node(node_id)
        # get info from node, and delete it from graph
        streaming_input, streaming_output, configuration_inputs = \
            node.get_single_streaming_input_single_output_and_configuration_inputs_of_node_for_parallelization()
        original_cmd_invocation_with_io_vars = node.cmd_invocation_with_io_vars

        prev_nodes = self.get_previous_nodes(node_id)
        first_pred_node, first_pred_cmd_inv = self.get_first_previous_node_and_first_previous_cmd_invocation(prev_nodes)

        # remove node to be parallelized
        self.remove_node(node_id) # remove it here already as as we need to remove edge end points ow. to avoid disconnecting graph to avoid disconnecting graph

        if len(prev_nodes) == 1 and isinstance(first_pred_node, r_merge.RMerge):
            # can be fused
            self.remove_node(prev_nodes[0]) # also sets respective edge to's and from's to None
            in_mapper_ids = first_pred_cmd_inv.operand_list
        else: # cannot be fused so introduce splitter
            # splitter
            round_robin_splitter_generator = lambda input_id, output_ids: r_split.make_r_split(input_id, output_ids, r_split_batch_size)
            out_split_ids = self.introduce_splitter(round_robin_splitter_generator, fan_out, fileIdGen, streaming_input)
            in_mapper_ids = out_split_ids

        # mappers
        out_mapper_ids = self.introduce_mappers(fan_out, fileIdGen, in_mapper_ids, original_cmd_invocation_with_io_vars,
                                                parallelizer)

        # aggregator
        self.introduce_aggregator_for_round_robin(out_mapper_ids, parallelizer, streaming_output)

    def apply_round_robin_with_unwrap_flag_parallelization_to_node(self, node_id, parallelizer, fileIdGen, fan_out,
                                                  batch_size, no_cat_split_vanish, r_split_batch_size):
        # round robin with unwrap flag is an inferred parallelizer which ensures that
        # the command is commutative and has an aggregator for consecutive chunks;
        # thus we can check whether we can re-open a previous "RR"-parallelization ending with `r_merge`
        node = self.get_node(node_id)
        streaming_input, streaming_output, configuration_inputs = \
            node.get_single_streaming_input_single_output_and_configuration_inputs_of_node_for_parallelization()
        original_cmd_invocation_with_io_vars = node.cmd_invocation_with_io_vars

        prev_nodes = self.get_previous_nodes(node_id)
        first_pred_node, first_pred_cmd_inv = self.get_first_previous_node_and_first_previous_cmd_invocation(prev_nodes)

        # remove node to be parallelized
        self.remove_node(node_id) # remove it here already as as we need to remove edge end points ow. to avoid disconnecting graph to avoid disconnecting graph

        if len(prev_nodes) == 1 and isinstance(first_pred_node, r_merge.RMerge):
                              # and node.is_commutative(): implied by how this kind of splitter is inferred
            self.remove_node(prev_nodes[0]) # also sets respective edge to's and from's to None

            in_unwrap_ids = first_pred_cmd_inv.operand_list
            out_unwrap_ids = self.introduce_unwraps(fileIdGen, in_unwrap_ids)
            in_mapper_ids = out_unwrap_ids
        else:
            # splitter
            round_robin_with_unwrap_flag_splitter_generator = lambda input_id, output_ids: r_split.make_r_split_with_unwrap_flag(input_id, output_ids, r_split_batch_size)
            out_split_ids = self.introduce_splitter(round_robin_with_unwrap_flag_splitter_generator, fan_out, fileIdGen, streaming_input)
            in_mapper_ids = out_split_ids

        # mappers
        out_mapper_ids = self.introduce_mappers(fan_out, fileIdGen, in_mapper_ids, original_cmd_invocation_with_io_vars,
                                                parallelizer)

        in_aggregator_ids = out_mapper_ids
        out_aggregator_id = streaming_output
        self.introduce_aggregators_for_consec_chunks(fileIdGen, in_aggregator_ids,
                                                     original_cmd_invocation_with_io_vars, out_aggregator_id, parallelizer,
                                                     streaming_output)

    def apply_consecutive_chunks_parallelization_to_node(self, node_id, parallelizer, fileIdGen, fan_out,
                                                               batch_size, no_cat_split_vanish, r_split_batch_size):
        # check whether we can fuse with previous node's parallelization:
        # we can do so if the previous node's parallelization is the same, and the aggregator is concatenation
        # Assumption: it suffices to check that the previous node is an aggregator node of type concatenate
        #  as this is unique for consecutive chunk parallelization (for now, this is true)
        node = self.get_node(node_id)
        streaming_input, streaming_output, configuration_inputs = \
            node.get_single_streaming_input_single_output_and_configuration_inputs_of_node_for_parallelization()
        original_cmd_invocation_with_io_vars = node.cmd_invocation_with_io_vars

        prev_nodes = self.get_previous_nodes(node_id)
        first_pred_node, first_pred_cmd_inv = self.get_first_previous_node_and_first_previous_cmd_invocation(prev_nodes)

        # remove node to be parallelized
        self.remove_node(node_id) # remove it here already as as we need to remove edge end points ow. to avoid disconnecting graph to avoid disconnecting graph

        # TODO: change first check to first_pred_node and not cmd_inv
        if len(prev_nodes) == 1 and first_pred_cmd_inv.is_aggregator_concatenate():
            # can be fused
            self.remove_node(prev_nodes[0]) # also sets respective edge to's and from's to None
            in_mapper_ids = first_pred_cmd_inv.operand_list
        else: # cannot be fused so introduce splitter
            # splitter
            consec_chunks_splitter_generator = lambda input_id, output_ids: pash_split.make_split_file(input_id,
                                                                                                       output_ids)
            out_split_ids = self.introduce_splitter(consec_chunks_splitter_generator, fan_out, fileIdGen,
                                                    streaming_input)
            in_mapper_ids = out_split_ids

        # mappers
        out_mapper_ids = self.introduce_mappers(fan_out, fileIdGen, in_mapper_ids, original_cmd_invocation_with_io_vars,
                                                parallelizer)

        # aggregators
        in_aggregator_ids = out_mapper_ids
        out_aggregator_id = streaming_output
        self.introduce_aggregators_for_consec_chunks(fileIdGen, in_aggregator_ids,
                                                     original_cmd_invocation_with_io_vars, out_aggregator_id, parallelizer,
                                                     streaming_output)

    def get_first_previous_node_and_first_previous_cmd_invocation(self, prev_nodes):
        assert (len(prev_nodes) > 0)
        # get info about first one but also ensure that it is the only one if we fuse
        first_pred_id = prev_nodes[0]
        first_pred_node = self.get_node(first_pred_id)
        first_pred_cmd_inv = first_pred_node.cmd_invocation_with_io_vars
        return first_pred_node, first_pred_cmd_inv

    def introduce_splitter(self, splitter_generator, fan_out, fileIdGen, streaming_input):
        out_split_ids = self.generate_ephemeral_edges(fileIdGen, fan_out)
        splitter = splitter_generator(streaming_input, out_split_ids)
        self.set_edge_to(streaming_input, splitter.get_id())
        for out_split_id in out_split_ids:
            self.set_edge_from(out_split_id, splitter.get_id())
        self.add_node(splitter)
        return out_split_ids

    def introduce_mappers(self, fan_out, fileIdGen, in_mapper_ids, original_cmd_invocation_with_io_vars, parallelizer):
        out_mapper_ids = self.generate_ephemeral_edges(fileIdGen, fan_out)
        zip_mapper_in_out_ids = zip(in_mapper_ids, out_mapper_ids)
        all_mappers = []
        for (in_id, out_id) in zip_mapper_in_out_ids:
            # BEGIN: these 4 lines could be refactored to be a function in graph such that
            # creating end point of edges and the creation of edges is not decoupled
            mapper_cmd_inv = parallelizer.get_actual_mapper(original_cmd_invocation_with_io_vars, in_id, out_id)
            mapper = DFGNode.make_simple_dfg_node_from_cmd_inv_with_io_vars(mapper_cmd_inv)
            self.set_edge_to(in_id, mapper.get_id())
            self.set_edge_from(out_id, mapper.get_id())
            # END
            splitter = parallelizer.get_splitter()
            if splitter.is_splitter_round_robin():
                mapper_r_wrapped = r_wrap.wrap_node(mapper, self.edges)
                self.set_edge_to(in_id, mapper_r_wrapped.get_id())
                self.set_edge_from(out_id, mapper_r_wrapped.get_id())
                mapper = mapper_r_wrapped
            all_mappers.append(mapper)
        for new_node in all_mappers:
            self.add_node(new_node)
        return out_mapper_ids

    def introduce_unwraps(self, fileIdGen, in_unwrap_ids):
        unwrap_to_commutative_mappers_ids = self.generate_ephemeral_edges(fileIdGen, len(in_unwrap_ids))
        in_out_unwrap_ids = zip(in_unwrap_ids, unwrap_to_commutative_mappers_ids)
        for in_unwrap, out_unwrap in in_out_unwrap_ids:
            unwrap = r_unwrap.make_unwrap_node([in_unwrap], out_unwrap)
            self.add_node(unwrap)
            self.set_edge_to(in_unwrap, unwrap.get_id())  # from are still (wrapped) mappers
            self.set_edge_from(out_unwrap, unwrap.get_id())  # to will be set to mappers of current node
        in_mapper_ids = unwrap_to_commutative_mappers_ids
        return in_mapper_ids

    def introduce_aggregators_for_consec_chunks(self, fileIdGen, in_aggregator_ids,
                                                original_cmd_invocation_with_io_vars, out_aggregator_id, parallelizer,
                                                streaming_output):
        aggregator_spec = parallelizer.get_aggregator_spec()
        if aggregator_spec.is_aggregator_spec_concatenate() or aggregator_spec.is_aggregator_spec_custom_n_ary():
            aggregator_cmd_inv = parallelizer.get_actual_aggregator(original_cmd_invocation_with_io_vars,
                                                                    in_aggregator_ids, out_aggregator_id)
            aggregator = DFGNode.make_simple_dfg_node_from_cmd_inv_with_io_vars(aggregator_cmd_inv)
            for in_aggregator_id in in_aggregator_ids:
                self.set_edge_to(in_aggregator_id, aggregator.get_id())
            self.set_edge_from(streaming_output, aggregator.get_id())
            all_aggregators = [aggregator]
            ## Add the merge commands in the graph
            for new_node in all_aggregators:
                self.add_node(new_node)
        elif aggregator_spec.is_aggregator_spec_custom_2_ary():
            # TODO: we simplify and assume that every mapper produces a single output for now
            map_in_aggregator_ids = [[id] for id in in_aggregator_ids]
            # TODO: turn node into cmd_invocation_with_io_vars since this is the only thing required in this function
            self.create_generic_aggregator_tree(original_cmd_invocation_with_io_vars, parallelizer, map_in_aggregator_ids, out_aggregator_id, fileIdGen)
        else:
            raise Exception("aggregator kind not yet implemented")

    def introduce_aggregator_for_round_robin(self, out_mapper_ids, parallelizer, streaming_output):
        aggregator_spec = parallelizer.get_aggregator_spec()
        if aggregator_spec.is_aggregator_spec_concatenate():
            in_aggregator_ids = out_mapper_ids
            out_aggregator_id = streaming_output
            aggregator = r_merge.make_r_merge_node(in_aggregator_ids, out_aggregator_id)
            for in_aggregator_id in in_aggregator_ids:
                self.set_edge_to(in_aggregator_id, aggregator.get_id())
            self.set_edge_from(streaming_output, aggregator.get_id())
            all_aggregators = [aggregator]
            ## Add the aggregator node(s) in the graph
            for new_node in all_aggregators:
                self.add_node(new_node)
        else:
            # TODO: this is where the other cases for aggregators need to be added
            pass


    ## This function parallelizes a merger followed by a parallelizable node
    ##
    ## There are several combinations that it can handle:
    ##   1. cat -> parallelizable node
    ##   2. r_merge -> stateless node without conf_input
    ##   3. r_merge -> commutative pure parallelizable node 
    ##
    ## 1. cat followed by a parallelizable node
    ##
    ##    (conf_input) ----+
    ##                      \
    ##    (in1) --- cat ---- node ---(out)---
    ##             /
    ##    (in2) --+
    ##
    ## is transformed to
    ##
    ##    (conf_input) -- tee ------------+
    ##                       \             \
    ##                        \   (in1) --- node --- agg ---(out)---
    ##                         \                    /
    ##                (in2) --- node --------------+
    ##
    ## where edges are named with parenthesis and nodes are named without them.
    ##
    ## 2. r_merge followed by a stateless node without conf_input
    ##
    ## TODO: Add visual representation
    ##
    ## In this case the stateless command is wrapped with wrap so we cannot actually tee the input (since we do not know apriori how many forks we have).
    ## However, we can actually write it to a file (not always worth performance wise) and then read it from all at once.
    ##
    ##
    ## TODO: Eventually delete the fileIdGen from here and always use the graph internal one.
    ##
    ## TODO: Eventually this should be tunable to not happen for all inputs (but maybe for less)
    def parallelize_node(self, node_id, fileIdGen):
        assert(False)
        node = self.get_node(node_id)
        # BEGIN ANNO
        # OLD
        # assert(node.is_parallelizable())
        # NEW
        rr_parallelizer_list = [parallelizer for parallelizer in node.parallelizer_list if parallelizer.splitter.is_splitter_round_robin()]
        assert(len(rr_parallelizer_list) == 1)
        rr_parallelizer = rr_parallelizer_list[0]
        # to have this info later when the merger is created in a reduce tree
        node.set_used_parallelizer(rr_parallelizer)
        # END ANNO

        ## Initialize the new_node list
        new_nodes = []

        ## Identify the previous merger node (cat or r_merge)
        ##
        ## TODO: This should also work for no cat (all inputs are part of the node)
        node_input_ids = node.get_standard_inputs()
        assert(len(node_input_ids) == 1)
        node_input_id = node_input_ids[0]
        previous_node_id = self.edges[node_input_id][1]
        previous_node = self.get_node(previous_node_id)
        assert(isinstance(previous_node, Cat)
               or isinstance(previous_node, r_merge.RMerge))

        ## Determine if the previous node is r_merge to determine which of the three parallelization cases to follow
        r_merge_flag = isinstance(previous_node, r_merge.RMerge)

        ## If the previous node of r_merge is an r_split, then we need to replace it with -r,
        ## instead of doing unwraps.
        if(r_merge_flag):
            assert(False)
            assert(isinstance(previous_node, r_merge.RMerge))
            r_merge_prev_node_ids = self.get_previous_nodes(previous_node_id)

            ## If all the previous nodes are r_split this means that they are the same
            ##
            ## Q: Could that ever not be true?
            ##
            ## TODO: If we ever want to measure the benefit from this optimization we need
            ##       to make a conjunction in this flag here.
            r_split_before_r_merge_opt_flag = all([isinstance(self.get_node(node_id), r_split.RSplit)
                                                   for node_id in r_merge_prev_node_ids])

            ## If r_split was right before the r_merge, and the node is pure parallelizable,
            ## this means that we will not add unwraps, and therefore we need to add the -r flag to r_split.
            if (r_split_before_r_merge_opt_flag
                and node.is_pure_parallelizable()):
                assert(node.is_commutative())
                r_split_id = r_merge_prev_node_ids[0]
                r_split_node = self.get_node(r_split_id)

                ## Add -r flag in r_split
                r_split_node.add_r_flag()
        else:
            r_split_before_r_merge_opt_flag = False


        ## Identify the parallel inputs, each of which will be given to a different copy of the node.
        parallel_input_ids = previous_node.get_input_list()
        parallelism = len(parallel_input_ids)

        ## Identify the output.
        node_output_edge_ids = node.get_output_list()
        assert(len(node_output_edge_ids) == 1)
        node_output_edge_id = node_output_edge_ids[0]

        ## Remove the original node and the cat node before it
        ## This also unplugs all the inputs
        self.remove_node(node_id)
        self.remove_node(previous_node_id)

        ## TODO: This does not work at the moment. There seem to be some issues with tee.
        ##       It probably has to do with a misunderstanding of how configuration inputs work
        ## Unplug the configuration inputs from the node and tee it
        parallel_configuration_ids = [[] for _ in range(parallelism)]
        node_conf_inputs = node.get_configuration_inputs()
        for conf_edge_id in node_conf_inputs:
            assert(False)
            ## TODO: For now this does not work for r_merge
            assert(not r_merge_flag)
            # self.set_edge_to(conf_edge_id, None)
            tee_id = self.tee_edge(conf_edge_id, parallelism, fileIdGen)
            tee_node = self.get_node(tee_id)
            for i in range(parallelism):
                # TODO outputs probably non-existent
                parallel_configuration_ids[i].append(tee_node.outputs[i])

        ## Create a temporary output edge for each parallel command.
        # BEGIN ANNO
        # OLD
        # map_output_fids = node.get_map_output_files(parallel_input_ids, fileIdGen)
        # NEW (added parameter)
        map_output_fids = get_map_output_files(node, parallel_input_ids, fileIdGen, rr_parallelizer)
        # END ANNO
        assert(len(map_output_fids) == len(parallel_input_ids))

        all_map_output_ids = []
        ## For each parallel input, create a parallel command
        for index in range(parallelism):
            ## Gather inputs and outputs
            conf_ins = parallel_configuration_ids[index]
            assert(len(conf_ins) == 0)
            standard_in = parallel_input_ids[index]
            new_inputs = (conf_ins, [standard_in])
            map_output_fid = map_output_fids[index]
            if(not isinstance(map_output_fid, list)):
                output_fid_list = [map_output_fid]
            else:
                output_fid_list = map_output_fid
            new_output_ids = [fid.get_ident() for fid in output_fid_list]
            all_map_output_ids.append(new_output_ids)

            ## Add the map_output_edges to the graph
            for output_fid in output_fid_list:
                self.add_edge(output_fid)

            ## If the previous merger is r_merge we need to put wrap around the nodes
            ## or unwrap before a commutative command
            if(r_merge_flag is True):
                assert(False)
                ## For stateless nodes we are in case (2) and we wrap them
                if (node.is_stateless()):
                    parallel_node = make_wrap_map_node(node, new_inputs, new_output_ids)
                    self.add_node(parallel_node)
                else:
                    ## If we have a pure parallelizable node, then we have to unwrap, before parallelizing the node.
                    ##
                    ## This can only work if the node is actually commutative
                    assert(node.is_pure_parallelizable())
                    assert(is_single_input(new_inputs))
                    assert(node.is_commutative())

                    ## Optimization: If the node before r_merge is an r_split, then we
                    ##               don't need to add unwrap, and we can just add -r to r_split.
                    if(r_split_before_r_merge_opt_flag):
                        parallel_node = make_map_node(node, new_inputs, new_output_ids)
                        self.add_node(parallel_node)
                    else:
                        ## Make the edge between unwrap and the command
                        unwrap_output_fid = fileIdGen.next_ephemeral_file_id()
                        unwrap_output_id = unwrap_output_fid.get_ident()
                        self.add_edge(unwrap_output_fid)

                        ## TODO: Make an unwrap node and create new inputs
                        unwrap_node = r_unwrap.make_unwrap_node(new_inputs, unwrap_output_id)
                        self.add_node(unwrap_node)
                        self.set_edge_from(unwrap_output_id, unwrap_node.get_id())

                        parallel_node_inputs = ([], [unwrap_output_id])
                        parallel_node = make_map_node(node, parallel_node_inputs, new_output_ids)
                        self.add_node(parallel_node)
                        self.set_edge_to(unwrap_output_id, parallel_node.get_id())

                        ## Note: unwrap needs to be set as the parallel node since below the inputs are set to point to it.
                        parallel_node = unwrap_node
            else:
                ## If we are working with a `cat` (and not an r_merge), then we just make a parallel node
                parallel_node = make_map_node(node, new_inputs, new_output_ids, rr_parallelizer)
                self.add_node(parallel_node)

            parallel_node_id = parallel_node.get_id()

            ## Set the to of all input edges
            for conf_in in conf_ins:
                assert(False)
                self.set_edge_to(conf_in, parallel_node_id)
            self.set_edge_to(standard_in, parallel_node_id)


        if (node.com_category == "stateless"):
            if(r_merge_flag is True):
                assert(False)
                new_merger = r_merge.make_r_merge_node(flatten_list(all_map_output_ids), node_output_edge_id)
            else:
                # BEGIN ANNO
                # OLD
                # new_merger = make_cat_node(flatten_list(all_map_output_ids), node_output_edge_id)
                # NEW
                new_merger = get_aggregator_as_dfg_node_from_node(node, rr_parallelizer, flatten_list(all_map_output_ids), [node_output_edge_id])
                # END ANNO

            self.add_node(new_merger)
            new_nodes.append(new_merger)
            self.set_edge_from(node_output_edge_id, new_merger.get_id())
        else:
            ## TODO: Create an aggregator here. At the moment it happens in `pash_runtime.py`.
            pass

        return new_nodes, all_map_output_ids

    ## Replicates an edge using tee and returns the new node_id.
    def tee_edge(self, edge_id, times, fileIdGen):
        ## Assert that the edge is unplugged
        assert(self.edges[edge_id][2] is None)

        output_fids = [fileIdGen.next_ephemeral_file_id() for _ in range(times)]
        output_ids = [fid.get_ident() for fid in output_fids]

        ## Create the tee node
        new_node = make_tee(edge_id, output_ids)
        new_node_id = new_node.get_id()

        ## Rewire the dfg
        for edge_fid in output_fids:
            self.add_from_edge(new_node_id, edge_fid)
        self.add_node(new_node)
        self.set_edge_to(edge_id, new_node_id)
        
        return new_node_id
        
    def generate_graphviz(self):
        ## TODO: It is unclear if importing in here (instead of in general)
        ##       improves startup cost of the pash_runtime when not using graphviz.
        import graphviz

        dot = graphviz.Digraph()

        ## TODO: Could use subgraph for distribution etc

        ## First generate all nodes
        for node_id, node in self.nodes.items():
            dot = node.add_dot_node(dot, node_id)

        ## (I/O) File nodes should be boxes
        dot.attr('node', shape='box')

        ## Then generate all edges and input+output files
        for fid, from_node, to_node in self.edges.values():
            ## This means that this file is an ending or starting file
            if from_node is None or to_node is None:
                ## Sometimes some source or sink files might be ephemeral
                ## TODO: We should investigate why this happens
                if fid.has_file_resource():
                    label = fid.serialize()
                    node_id = f'file-{str(fid.get_ident())}'
                    dot.node(node_id, label)

                    if from_node is None:
                        dot.edge(node_id, str(to_node))
                    else:
                        dot.edge(str(from_node), node_id)
            else:
                dot.edge(str(from_node), str(to_node))

        return dot

    ## TODO: Also it should check that there are no unreachable edges

    def edge_node_consistency(self):
        ## Check if edges and nodes are consistent
        for edge_id, (_, from_node_id, to_node_id) in self.edges.items():
            if (not from_node_id is None):
                from_node = self.get_node(from_node_id)
                if(not (edge_id in from_node.get_output_list())):
                    log("Consistency Error: Edge id:", edge_id, "is not in the node outputs:", from_node)
                    return False
            if (not to_node_id is None):
                to_node = self.get_node(to_node_id)
                if(not (edge_id in to_node.get_input_list())):
                    log("Consistency Error: Edge id:", edge_id, "is not in the node inputs:", to_node)
                    return False

        for node_id, node in self.nodes.items():
            for edge_id in node.get_input_list():
                _, _, to_node_id = self.edges[edge_id]
                if(not (to_node_id == node_id)):
                    log("Consistency Error: The to_node_id of the input_edge:", edge_id, "of the node:", node, "is equal to:", to_node_id)
                    return False
            for edge_id in node.get_output_list():
                _, from_node_id, _ = self.edges[edge_id]
                if(not (from_node_id == node_id)):
                    log("Consistency Error: The from_node_id of the output_edge:", edge_id, "of the node:", node, "is equal to:", from_node_id)
                    return False

        return True

    ## This function checks whether an IR is valid -- that is, if it
    ## has at least one node, and stdin, stdout set to some non-null
    ## file identifiers.
    def valid(self):
        return (len(self.nodes) > 0 and
                self.edge_node_consistency() and
                (not self.is_in_background()
                  or (self.get_stdin() is None)))
                ## The following is not true. Background IRs should not have stdin, but they can have stdout.
                #   and self.get_stdout() is None)))
                ## The following is not true. A DFG might not have an stdin
                #  or (not self.is_in_background()
                #      and not self.get_stdin() is None 
                #      and not self.get_stdout() is None)))

    ## This is a function that creates a reduce tree for a given node
    def create_generic_aggregator_tree(self, cmd_invocation_with_io_vars, parallelizer, input_ids_for_aggregators, out_aggregator_id, fileIdGen):
        def function_to_get_binary_aggregator(in_ids, out_ids):
            assert(len(out_ids) == 1)
            aggregator_cmd_inv = parallelizer.get_actual_aggregator(cmd_invocation_with_io_vars, in_ids, out_ids[0])
            aggregator = DFGNode.make_simple_dfg_node_from_cmd_inv_with_io_vars(aggregator_cmd_inv)
            return aggregator
        ## The Aggregator node takes a sequence of input ids and an output id
        all_aggregators, new_edges, final_output_id = self.create_reduce_tree(lambda in_ids, out_ids: function_to_get_binary_aggregator(in_ids, out_ids),
                                    input_ids_for_aggregators, fileIdGen)
        ## Add the edges in the graph
        self.add_edges(new_edges)
        ## Add the merge commands in the graph
        for new_node in all_aggregators:
            self.add_node(new_node)

        ## Replace the previous final_output_id with the previous id
        node_output_edge_id = out_aggregator_id
        final_merge_node_id = self.edges[final_output_id][1]
        final_merge_node = self.get_node(final_merge_node_id)
        final_merge_node.replace_edge(final_output_id, node_output_edge_id)
        self.set_edge_from(node_output_edge_id, final_merge_node_id)
        self.set_edge_from(final_output_id, None)

    ## This function creates the reduce tree. Both input and output file
    ## ids must be lists of lists, as the input file ids and the output
    ## file ids might contain auxiliary files.
    def create_reduce_tree(self, init_func, input_ids, fileIdGen):
        tree = []
        new_edges = []
        curr_ids = input_ids
        while(len(curr_ids) > 1):
            new_level, curr_ids, new_fids = self.create_reduce_tree_level(init_func, curr_ids, fileIdGen)
            tree += new_level
            new_edges += new_fids

        # Find the final output    (provided with parameter)
        final_output_id = curr_ids[0][0]

        ## Drain the final auxiliary outputs
        final_auxiliary_outputs = curr_ids[0][1:]
        drain_fids = [fileIdGen.next_file_id()
                      for final_auxiliary_output in final_auxiliary_outputs]
        for drain_fid in drain_fids:
            drain_fid.set_resource(FileResource(Arg(string_to_argument('/dev/null'))))
            new_edges.append(drain_fid)
        drain_ids = [fid.get_ident() for fid in drain_fids]

        drain_cat_commands = [make_cat_node([final_auxiliary_output], drain_id)
                              for final_auxiliary_output, drain_id in zip(final_auxiliary_outputs, drain_ids)]
        return (tree + drain_cat_commands), new_edges, final_output_id

    @staticmethod
    ## This function creates a level of the reduce tree. Both input and
    ## output file ids must be lists of lists, as the input file ids and
    ## the output file ids might contain auxiliary files.
    def create_reduce_tree_level(init_func, input_ids, fileIdGen):
        if(len(input_ids) % 2 == 0):
            output_ids = []
            even_input_ids = input_ids
        else:
            output_ids = [input_ids[0]]
            even_input_ids = input_ids[1:]

        new_fids = []
        level = []
        for i in range(0, len(even_input_ids), 2):
            new_out_fids = [fileIdGen.next_ephemeral_file_id() for _ in input_ids[i]]
            new_fids += new_out_fids
            new_out_ids = [fid.get_ident() for fid in new_out_fids]
            output_ids.append(new_out_ids)
            new_node = IR.create_reduce_node(init_func, even_input_ids[i:i+2], new_out_ids)
            level.append(new_node)
        return (level, output_ids, new_fids)

    @staticmethod
    ## This function creates one node of the reduce tree
    def create_reduce_node(init_func, input_ids, output_ids):
        return init_func(flatten_list(input_ids), output_ids)
    # TODO: this is where we need to use our aggregator spec/node

