import copy
import json
import yaml
import os

from definitions.ir.arg import *
from definitions.ir.node import *
from definitions.ir.command import *
from definitions.ir.resource import *
from definitions.ir.nodes.cat import *

from command_categories import *
from union_find import *
from ir_utils import *
from util import *

import config

## Creates a file id for a given resource
def create_file_id_for_resource(resource, fileIdGen):
    file_id = create_split_file_id(resource.get_length(), fileIdGen)
    file_id.set_resource(resource)
    return file_id

## Creates a file id that has a given maximum length
def create_split_file_id(batch_size, fileIdGen):
    file_id = fileIdGen.next_file_id()
    file_id.set_max_length(batch_size)
    return file_id

class FileIdGen:
    def __init__(self, next = 0, prefix = ""):
        self.next = next + 1
        self.prefix = prefix

    def next_file_id(self):
        fileId = FileId(self.next, self.prefix)
        self.next += 1
        return fileId

def create_command_assign_file_identifiers(ast, fileIdGen, command, options,
                                           stdin=None, stdout=None, redirections=[]):
    in_stream, out_stream, opt_indices = find_command_input_output(command, options, stdin, stdout)
    category = find_command_category(command, options)
    ## TODO: Maybe compile command instead of just formatting it.
    if(format_arg_chars(command) == "cat"):
        command = Cat(ast, command, options, in_stream, out_stream,
                      opt_indices, category, stdin, stdout, redirections)
    else:
        command = Command(ast, command, options, in_stream, out_stream,
                          opt_indices, category, stdin, stdout, redirections)

    ## The options that are part of the input and output streams must
    ## be swapped with file identifiers. This means that each file
    ## identifier must have a unique resource that it points to.
    for opt_or_ch in in_stream:
        new_fid = replace_file_arg_with_id(opt_or_ch, command, fileIdGen)
        command.set_file_id(opt_or_ch, new_fid)

    for opt_or_ch in out_stream:
        new_fid = replace_file_arg_with_id(opt_or_ch, command, fileIdGen)
        command.set_file_id(opt_or_ch, new_fid)

    return command

def replace_file_arg_with_id(opt_or_channel, command, fileIdGen):
    fid_or_resource = command.get_file_id(opt_or_channel)
    ## If the file is not a FileId, then it is some argument. We
    ## create a file identifier, and replace it with that, and
    ## make sure that the file identifier points to the argument.
    if (not isinstance(fid_or_resource, FileId)):
        return create_file_id_for_resource(Resource(fid_or_resource), fileIdGen)
    else:
        return fid_or_resource


def make_split_files(in_fid, fan_out, batch_size, fileIdGen):
    assert(fan_out > 1)
    ## Generate the split file ids
    out_fids = [fileIdGen.next_file_id() for i in range(fan_out)]
    split_com = make_split_file(in_fid, out_fids, batch_size)
    return [split_com], out_fids

## TODO: Make a proper splitter subclass of Node
def make_split_file(in_fid, out_fids, batch_size):
    ## TODO: I probably have to give the file names as options to the command to.
    options = [string_to_argument(str(batch_size))] + out_fids
    opt_indices = [("option", 0)]
    stdout_indices = [("option", i+1) for i in range(len(out_fids))]
    command = Command(None, # TODO: Make a proper AST
                      string_to_argument("split_file"),
                      options,
                      ["stdin"],
                      stdout_indices,
                      opt_indices,
                      None, # TODO: Category?
                      in_fid)
    return command

## This function gets a file identifier and returns the maximum among
## its, and its parents identifier (parent regarding Union Find)
def get_larger_file_id_ident(file_ids):
    return max([max(fid.get_ident(), Find(fid).get_ident())
                for fid in file_ids])

## Note: This might need more information. E.g. all the file
## descriptors of the IR, and in general any other local information
## that might be relevant.
class IR:

    ## TODO: Improve the rerpesentation and keep the input and output fids in
    ##       some structure.

    ## IR Assumptions:
    ##
    ## - Each node has a list of incoming files in order of
    ##   consumption.
    ##
    ## - If two nodes have the same file as output, then they both
    ##   write to it concurrently.
    def __init__(self, nodes, stdin = [], stdout = [], background = False):
        self.nodes = nodes
        self.stdin = stdin
        self.stdout = stdout
        self.background = background

    def __repr__(self):
        output = "(|-{} IR: {} {}-|)".format(self.stdin, self.nodes, self.stdout)
        return output

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

        ## Redirect inputs
        in_fids = self.all_input_fids()
        ## TODO: Remove this assert, extending input to more than stdin
        assert(len(in_fids) <= 1)
        for in_fid in in_fids:
            file_to_redirect_to = in_fid.to_ast()
            redirect_stdin_script = os.path.join(config.PASH_TOP, config.config['runtime']['redirect_stdin_binary'])
            com_args = [string_to_argument('source'), string_to_argument(redirect_stdin_script), file_to_redirect_to]
            com = make_command(com_args)
            asts.append(com)

        ## Make the dataflow graph
        for node in self.nodes:
            node_ast = node.to_ast(drain_streams)
            asts.append(make_background(node_ast))
        
        ## Redirect outputs
        out_fids = self.all_output_fids()
        ## TODO: Make this work for more than one output. 
        ##       For now it is fine to only have stdout as output
        ##
        ## TODO: Make this not cat if the output is a real file.
        assert(len(out_fids) == 1)
        com_args = [string_to_argument('cat'), out_fids[0].to_ast()]
        node = make_background(make_command(com_args))
        asts.append(node)

        return asts


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
        self.nodes += other.nodes

        ## This combines the two IRs by adding all of the nodes
        ## together, and by union-ing the stdout of the first with the
        ## stdin of the second.
        ##
        ## Question: What happens if one of them is NULL. This
        ##           shouldn't be the case after we check that
        ##           both self and other are not empty.
        assert(len(self.stdout) == 1)
        assert(len(self.stdin) == 1)
        self.stdout[0].union(other.stdin[0])
        self.stdout = other.stdout

        ## Note: The ast is not extensible, and thus should be
        ## invalidated if an operation happens on the IR
        self.ast = None

    def union(self, other):
        assert(self.valid())
        assert(other.valid())
        assert(self.is_in_background())
        self.nodes += other.nodes

        ## This combines two IRs where at least the first one is in
        ## background. This means that the stdin, stdout are those of
        ## the second (or None if both are in background). Also if
        ## both are in background, their union is also in background.

        ## TODO: Handle any redirections

        ## If one of them is not in the background, then the whole
        ## thing isn't.
        if (not other.is_in_background()):
            self.set_background(other.is_in_background())
            self.stdin = other.stdin
            self.stdout = other.stdout

        ## Note: The ast is not extensible, and thus should be
        ## invalidated if an operation happens on the IR
        self.ast = None

        ## TODO: Handle connections of common files (pipes, etc)
        self.combine_common_files()


    ## Combines (unions) files that refer to the same resource.
    ##
    ## WARNING: This assumes that comparing file names statically
    ## (syntactically) for equality, implies semantic
    ## equality. However this is not true. There are cases where
    ## different identifiers could refer to the same file.
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

        ## TODO: (Maybe) Signal an error if an output is an output in
        ## more than one node

        ## For all inputs of all nodes, check if they are the output
        ## of exactly one other node.
        for node in self.nodes:
            in_stream_with_resources = [file_in for file_in in node.get_input_file_ids()
                                        if file_in.has_resource()]
            for file_in in in_stream_with_resources:
                in_resource = file_in.get_resource()
                number_of_out_resources = 0
                for node2 in self.nodes:
                    out_stream_with_resources = [file_out for file_out in node.get_output_file_ids()
                                                 if file_out.has_resource()]
                    for file_out in out_stream_with_resources:
                        if (in_resource == out_resource):
                            number_of_out_resources += 1
                            file_in.union(file_out)
                assert(number_of_out_resources <= 1)

    ## Returns all the file identifiers in the IR.
    def all_fids(self):
        all_file_ids = []
        for node in self.nodes:
            ## Gather all fids that this node uses.
            input_pipes = [fid for fid in node.get_input_file_ids()]
            output_pipes = [fid for fid in node.get_output_file_ids()]
            all_file_ids += input_pipes + output_pipes

        ## Remove duplicates
        ## TODO: This should normally happen without comparing strings, 
        ##       but there seems to be some problem now and some fifo is 
        ##       considered both an argument and an fid.
        unique_file_ids = []
        seen_file_id_strings = set()
        for fid in all_file_ids:
            fid_string = fid.serialize()
            if(not fid_string in seen_file_id_strings):
                unique_file_ids.append(fid)
                seen_file_id_strings.add(fid_string)
        # all_file_ids = list(set(all_file_ids))
        return unique_file_ids

    ## Returns all input fids of the IR
    def all_input_fids(self):
        flat_stdin = [Find(fid) for fid in self.stdin]
        return flat_stdin

    ## Returns all output fids of the IR
    def all_output_fids(self):
        flat_stdout = [Find(fid) for fid in self.stdout]
        return flat_stdout

    ## Returns the sources of the IR (i.e. the nodes that has no
    ## incoming edge)
    def source_nodes(self):
        sources = [node for node in self.nodes if not self.has_incoming_edge(node)]
        return sources

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

    def get_next_nodes_and_edges(self, node):
        next_nodes_and_edges = []
        for outgoing_fid in node.get_output_file_ids():
            for other_node in self.nodes:
                ## Note: What if other_node == node?
                if (not outgoing_fid.find_fid_list(other_node.get_input_file_ids()) is None):
                    next_nodes_and_edges.append((other_node, outgoing_fid))
        return next_nodes_and_edges

    def get_next_nodes(self, node):
        return [node for node, edge in self.get_next_nodes_and_edges(node)]

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
        max_id = 0
        max_id = max(get_larger_file_id_ident(self.stdin), max_id)
        max_id = max(get_larger_file_id_ident(self.stdout), max_id)
        for node in self.nodes:
            node_file_ids = node.get_input_file_ids() + node.get_output_file_ids()
            for file_id in node_file_ids:
                max_id = max(get_larger_file_id_ident([file_id]), max_id)
        return FileIdGen(max_id)



    def remove_node(self, node):
        self.nodes.remove(node)

    def add_node(self, node):
        self.nodes.append(node)

    ## Note: We assume that the lack of nodes is an adequate condition
    ##       to check emptiness.
    def empty(self):
        return (len(self.nodes) == 0)

    ## This function checks whether an IR is valid -- that is, if it
    ## has at least one node, and stdin, stdout set to some non-null
    ## file identifiers.
    def valid(self):
        return (len(self.nodes) > 0 and
                ((self.is_in_background()
                  and len(self.stdin) == 0
                  and len(self.stdout) == 0)
                 or (not self.is_in_background()
                     and len(self.stdin) > 0
                     and len(self.stdout) > 0)))

