import copy
import json
import yaml
import os

from definitions.ir.arg import *
from definitions.ir.node import *
from definitions.ir.command import *
from definitions.ir.resource import *
from definitions.ir.nodes.cat import *
from definitions.ir.nodes.pash_split import *

from command_categories import *
from union_find import *
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

## This function creates a DFG with a single node given a command.
def compile_command_to_DFG(fileIdGen, command, options,
                           redirections=[]):
    ## TODO: There is no need for this redirection here. We can just straight 
    ##       come up with inputs, outputs, options
    in_stream, out_stream, opt_indices = find_command_input_output(command, options)
    # log("Opt indices:", opt_indices, "options:", options)
    category = find_command_category(command, options)

    ## TODO: Make an empty IR and add edges and nodes incrementally (using the methods defined in IR).

    dfg_edges = {}
    ## Add all inputs and outputs to the DFG edges
    dfg_inputs = []
    for opt_or_fd in in_stream:
        fid = get_option_or_fd(opt_or_fd, options, fileIdGen)
        fid_id = fid.get_ident()
        dfg_edges[fid_id] = (fid, None, None)
        dfg_inputs.append(fid_id)

    dfg_outputs = []
    for opt_or_fd in out_stream:
        fid = get_option_or_fd(opt_or_fd, options, fileIdGen)
        fid_id = fid.get_ident()
        dfg_edges[fid_id] = (fid, None, None)
        dfg_outputs.append(fid_id)

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
                       ## TODO: We don't really need to pass category for Cat
                       com_category,
                       com_options=dfg_options,
                       com_redirs=com_redirs,
                       com_assignments=com_assignments)
        ## TODO: Make Cat a subtype of DFGNode
        # command = Cat(ast, command, options, in_stream, out_stream,
        #               opt_indices, category, stdin, stdout, redirections)
    else:
        ## Assume: Everything must be completely expanded
        ## TODO: Add an assertion about that.
        dfg_node = DFGNode(dfg_inputs, 
                           dfg_outputs, 
                           com_name,
                           com_category,
                           com_options=dfg_options,
                           com_redirs=com_redirs,
                           com_assignments=com_assignments)
    
    if(not dfg_node.is_at_most_pure()):
        raise ValueError()

    node_id = id(dfg_node)

    ## Assign the from, to node in edges
    for fid_id in dfg_inputs:
        fid, from_node, to_node = dfg_edges[fid_id]
        assert(to_node is None)
        dfg_edges[fid_id] = (fid, from_node, node_id)
    
    for fid_id in dfg_outputs:
        fid, from_node, to_node = dfg_edges[fid_id]
        assert(from_node is None)
        dfg_edges[fid_id] = (fid, node_id, to_node)
    
    dfg_nodes = {node_id : dfg_node}
    dfg = IR(dfg_nodes, dfg_edges)
    return dfg



## TODO: Delete all unneccessary functions.

# def replace_file_arg_with_id(opt_or_channel, command, fileIdGen):
#     fid_or_resource = command.get_file_id(opt_or_channel)
#     ## If the file is not a FileId, then it is some argument. We
#     ## create a file identifier, and replace it with that, and
#     ## make sure that the file identifier points to the argument.
#     if (not isinstance(fid_or_resource, FileId)):
#         return create_file_id_for_resource(Resource(fid_or_resource), fileIdGen)
#     else:
#         return fid_or_resource


def make_split_files(in_fid, fan_out, batch_size, fileIdGen):
    assert(fan_out > 1)
    ## Generate the split file ids
    out_fids = [fileIdGen.next_file_id() for i in range(fan_out)]
    split_com = make_split_file(in_fid, out_fids, batch_size)
    return [split_com], out_fids

## This function gets a file identifier and returns the maximum among
## its, and its parents identifier (parent regarding Union Find)
def get_larger_file_id_ident(file_ids):
    return max([max(fid.get_ident(), Find(fid).get_ident())
                for fid in file_ids])


## TODO: Revise comments and all methods
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
        for _node_id, node in self.nodes.items():
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
        assert(self.valid())
        assert(other.valid())

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
        for _node_id1, node1 in self.nodes.items():
            inputs_with_file_resource = [(id, fid) for id, fid in node1.get_input_ids_fids(self.edges)
                                         if fid.has_file_resource()]
            for id_in, fid_in in inputs_with_file_resource:
                in_resource = fid_in.get_resource()
                number_of_out_resources = 0
                for node_id2, node2 in self.nodes.items():
                    outputs_with_file_resource = [(id, fid) for id, fid in node2.get_output_ids_fids(self.edges)
                                                  if fid.has_file_resource()]
                    for id_out, fid_out in outputs_with_file_resource:
                        out_resource = fid_out.get_resource()
                        if (in_resource == out_resource):
                            number_of_out_resources += 1
                            ## They point to the same File resource so we need to unify their fids
                            self.nodes[node_id2].replace_edge(id_out, id_in)

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

    ## TODO: Delete this
    ## This function returns whether a node has an incoming edge in an IR
    ##
    ## WARNING: At the moment is is extremely naive and slow.
    def has_incoming_edge(self, node):
        for incoming_fid in node.get_input_file_ids():
            for other_node in self.nodes:
                ## Note: What if other_node == node?
                if (not incoming_fid.find_fid_list(other_node.get_output_file_ids())
                    is None):
                    return True
        return False

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
    
    def get_node_input_ids_fids(self, node_id):
        node = self.get_node(node_id)
        return [(input_edge_id, self.edges[input_edge_id][0]) for input_edge_id in node.inputs]

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

    ## TODO: Delete this
    def get_previous_nodes_and_edges(self, node):
        previous_nodes_and_edges = []
        for incoming_fid in node.get_input_file_ids():
            previous_node = self.get_fids_previous_node_and_edge(incoming_fid)
            if (not previous_node is None):
                previous_nodes_and_edges.append((previous_node, incoming_fid))
        return previous_nodes_and_edges

    def get_fids_previous_node_and_edge(self, fid):
        for other_node in self.nodes:
            ## Note: What if other_node == node?
            if (not fid.find_fid_list(other_node.get_output_file_ids()) is None):
                return other_node
        return None

    def get_previous_nodes(self, node):
        return [node for node, edge in self.get_previous_nodes_and_edges(node)]

    ## This command gets all file identifiers of the graph, and
    ## returns a fileId generator that won't clash with the existing
    ## ones.
    def get_file_id_gen(self):
        max_id = max(self.edges.keys())
        return FileIdGen(max_id)

    def remove_node(self, node_id):
        node = self.nodes.pop(node_id)
        ## Remove the node in the edges dictionary
        for in_id in node.inputs:
            self.set_edge_to(in_id, None)
        
        for out_id in node.outputs:
            self.set_edge_from(out_id, None)


    def add_node(self, node):
        node_id = id(node)
        self.nodes[node_id] = node
        ## Add the node in the edges dictionary
        for in_id in node.inputs:
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

    ## TODO: Add a part of valid that checks if all nodes and edges are consistent.
    ## TODO: Also it should check that there are no unreachable edges

    ## This function checks whether an IR is valid -- that is, if it
    ## has at least one node, and stdin, stdout set to some non-null
    ## file identifiers.
    def valid(self):
        return (len(self.nodes) > 0 and
                (not self.is_in_background()
                  or (self.get_stdin() is None)))
                    ## The following is not true. Background IRs should not have stdin, but they can have stdout.
                    #   and self.get_stdout() is None)))
                 ## The following is not true. A DFG might not have an stdin
                #  or (not self.is_in_background()
                #      and not self.get_stdin() is None 
                #      and not self.get_stdout() is None)))

