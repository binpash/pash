import json
import yaml
import os

from definitions.ir.arg import *
from definitions.ir.dfg_node import *
from definitions.ir.file_id import *
from definitions.ir.resource import *
from definitions.ir.nodes.cat import *
from definitions.ir.nodes.bigram_g_map import *

import definitions.ir.nodes.pash_split as pash_split
import definitions.ir.nodes.r_merge as r_merge
import definitions.ir.nodes.r_split as r_split
import definitions.ir.nodes.r_wrap as r_wrap
import definitions.ir.nodes.r_unwrap as r_unwrap

from command_categories import *
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
        self.prefix = prefix

    def next_file_id(self):
        fileId = FileId(self.next, self.prefix)
        self.next += 1
        return fileId

    def next_ephemeral_file_id(self):
        fileId = FileId(self.next, self.prefix)
        self.next += 1
        fileId.make_ephemeral()
        return fileId

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

def find_input_edges(inputs, dfg_edges, options, fileIdGen):
    if(isinstance(inputs, list)):
        return create_edges_from_opt_or_fd_list(inputs, dfg_edges, options, fileIdGen)
    elif(isinstance(inputs, tuple)):
        config_inputs = create_edges_from_opt_or_fd_list(inputs[0], dfg_edges, options, fileIdGen)
        standard_inputs = create_edges_from_opt_or_fd_list(inputs[1], dfg_edges, options, fileIdGen)
        return (config_inputs, standard_inputs)
    else:
        raise NotImplementedError()

## This function creates a DFG with a single node given a command.
def compile_command_to_DFG(fileIdGen, command, options,
                           redirections=[]):
    ## TODO: There is no need for this redirection here. We can just straight 
    ##       come up with inputs, outputs, options
    inputs, out_stream, opt_indices = find_command_input_output(command, options)
    # log("Opt indices:", opt_indices, "options:", options)
    category = find_command_category(command, options)
    com_properties = find_command_properties(command, options)
    com_aggregator = find_command_aggregator(command, options)

    ## TODO: Make an empty IR and add edges and nodes incrementally (using the methods defined in IR).

    dfg_edges = {}
    ## Add all inputs and outputs to the DFG edges
    dfg_inputs = find_input_edges(inputs, dfg_edges, options, fileIdGen)
    dfg_outputs = create_edges_from_opt_or_fd_list(out_stream, dfg_edges, options, fileIdGen)

    com_name = Arg(command)
    com_category = category

    ## Get the options
    dfg_options = [get_option(opt_or_fd, options, fileIdGen)
                   for opt_or_fd in opt_indices]
    com_redirs = redirections
    ## TODO: Add assignments
    com_assignments = []

    ## TODO: Combine them both in a constructor that decided whether to instantiate Cat or DFGNode
    if(str(com_name) == "cat"):
        dfg_node = Cat(dfg_inputs,
                       dfg_outputs,
                       com_name,
                       ## TODO: We don't really need to pass category, name, or input_consumption for Cat
                       com_category,
                       com_options=dfg_options,
                       com_redirs=com_redirs,
                       com_assignments=com_assignments)
    else:
        ## Assume: Everything must be completely expanded
        ## TODO: Add an assertion about that.
        dfg_node = DFGNode(dfg_inputs, 
                           dfg_outputs, 
                           com_name,
                           com_category,
                           com_properties=com_properties,
                           com_aggregator=com_aggregator,
                           com_options=dfg_options,
                           com_redirs=com_redirs,
                           com_assignments=com_assignments)
    
    if(not dfg_node.is_at_most_pure()):
        raise ValueError()

    node_id = dfg_node.get_id()

    ## Assign the from, to node in edges
    for fid_id in dfg_node.get_input_list():
        fid, from_node, to_node = dfg_edges[fid_id]
        assert(to_node is None)
        dfg_edges[fid_id] = (fid, from_node, node_id)
    
    for fid_id in dfg_node.outputs:
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

## TODO: Move it somewhere else, but where?
def make_map_node(node, new_inputs, new_outputs):
    ## Some nodes have special map commands
    ##
    ## TODO: Make this more general instead of hardcoded
    if(str(node.com_name) == "bigrams_aux"):
        ## Ensure that the inputs have the correct size for this
        assert(len(new_inputs[0]) == 0)
        assert(len(new_inputs[1]) == 1)
        new_node = BigramGMap(new_inputs[1][0], new_outputs)
    else:
        new_node = node.copy()
        new_node.inputs = new_inputs
        new_node.outputs = new_outputs
    return new_node

## Makes a wrap node that encloses a map parallel node.
##
## At the moment it only works with one input and one output since wrap cannot redirect input in the command.
def make_wrap_map_node(node, new_inputs, new_outputs):
    # log("Inputs:", new_inputs)
    # log("Outputs:", new_outputs)
    assert(is_single_input(new_inputs))
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
        # log("Nodes:", self.nodes)
        # log("Edges:", self.edges)

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
                assert(stdout_id is None)
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

        ## Put the output node in the end for wait to work.
        for node_id in sink_node_ids:
            node = self.get_node(node_id)
            node_ast = node.to_ast(self.edges, drain_streams)
            asts.append(make_background(node_ast))

        return asts

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

    ## Returns the sources of the IR (i.e. the nodes that has no
    ## incoming edge)
    def source_nodes(self):
        sources = set()
        for _edge_fid, from_node, to_node in self.edges.values():
            if(from_node is None and not to_node is None):
                sources.add(to_node)
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
        output_edge_ids = self.nodes[node_id].outputs
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
        return [(output_edge_id, self.edges[output_edge_id][0]) for output_edge_id in node.outputs]

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
        
        for out_id in node.outputs:
            self.set_edge_from(out_id, None)


    def add_node(self, node):
        node_id = node.get_id()
        self.nodes[node_id] = node
        ## Add the node in the edges dictionary
        for in_id in node.get_input_list():
            self.set_edge_to(in_id, node_id)
        
        for out_id in node.outputs:
            self.set_edge_from(out_id, node_id)


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
        node = self.get_node(node_id)
        assert(node.is_parallelizable())

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
        node_output_edge_ids = node.outputs
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
            ## TODO: For now this does not work for r_merge
            assert(not r_merge_flag)
            # self.set_edge_to(conf_edge_id, None)
            tee_id = self.tee_edge(conf_edge_id, parallelism, fileIdGen)
            tee_node = self.get_node(tee_id)
            for i in range(parallelism):
                parallel_configuration_ids[i].append(tee_node.outputs[i])
        
        ## Create a temporary output edge for each parallel command.
        map_output_fids = node.get_map_output_files(parallel_input_ids, fileIdGen)

        all_map_output_ids = []
        ## For each parallel input, create a parallel command
        for index in range(parallelism):
            ## Gather inputs and outputs
            conf_ins = parallel_configuration_ids[index]
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
                parallel_node = make_map_node(node, new_inputs, new_output_ids)
                self.add_node(parallel_node)

            parallel_node_id = parallel_node.get_id()

            ## Set the to of all input edges
            for conf_in in conf_ins:
                self.set_edge_to(conf_in, parallel_node_id)
            self.set_edge_to(standard_in, parallel_node_id)


        if (node.com_category == "stateless"):
            if(r_merge_flag is True):
                new_merger = r_merge.make_r_merge_node(flatten_list(all_map_output_ids), node_output_edge_id)
            else:
                new_merger = make_cat_node(flatten_list(all_map_output_ids), node_output_edge_id)
            
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
        


    ## TODO: Also it should check that there are no unreachable edges

    def edge_node_consistency(self):
        ## Check if edges and nodes are consistent
        for edge_id, (_, from_node_id, to_node_id) in self.edges.items():
            if (not from_node_id is None):
                from_node = self.get_node(from_node_id)
                if(not (edge_id in from_node.outputs)):
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
            for edge_id in node.outputs:
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

