import copy
import json
import yaml
import os

from definitions.ir.arg import *
from definitions.ir.node import *
from definitions.ir.command import *
from definitions.ir.resource import *
from definitions.ir.nodes.cat import *

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
                      opt_indices, category, stdin, stdout)
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


## TODO: Make dictionary that holds command information (category,
## inputs-outputs, etc...)

## This function returns the input and output streams of a command.
##
## The input and output lists, contain tuples that refer to options:
## e.g. ("option", 0) or "stdin", "stdout" when they refer to stdin or
## stdout.
##
## At the moment it has just hardcoded knowledge of the inputs and
## outputs of several commands.
##
## By default they are the stdin and the stdout of the node, and they
## are only filled in for commands that we (or the developer) has
## specified a list of input resources that also contains files in the
## arguments.
def find_command_input_output(command, options, stdin, stdout):
    command_string = format_arg_chars(command)
    # print("Command to categorize:", command_string)

    assert(isinstance(command_string, str))

    ## TODO: Make a proper search that returns the command outputs and
    ## inputs. This is hardcoded and wrong
    print(" -- Warning: Argument inputs and outputs for: {} are hardcoded and possibly wrong"
          .format(command_string))

    if (command_string == "cat"):
        input_stream = [("option", i) for i in range(len(options))]
        return (input_stream, ["stdout"], [])
    elif (command_string == "comm"):
        return comm_input_output(options, stdin, stdout)
    elif (command_string == "diff"):
        return diff_input_output(options, stdin, stdout)
    else:
        opt_indices = [("option", i) for i in range(len(options))]
        return (["stdin"], ["stdout"], opt_indices)


## This functions finds and returns a string representing the command category
def find_command_category(command, options):
    command_string = format_arg_chars(command)
    print("Command to categorize:", command_string)

    assert(isinstance(command_string, str))

    ## TODO: Make a proper search that returns the command category
    print(" -- Warning: Category for: {} is hardcoded and possibly wrong".format(command_string))

    # NOTE order of class declaration in definition file is important, as it
    # dictates class precedence in the following search
    for command_class, commands in config.command_classes.items():
        command_list = list(map(get_command_from_definition, commands))

        if (command_string in command_list
            or command_string.split("/")[-1] in command_list):
            return command_class

    if command_string == 'comm':
        return is_comm_pure(options)

    return 'none'

## TODO: This is clearly incomplete
def is_comm_pure(options):
    first_opt = format_arg_chars(options[0])
    if(first_opt == "-13" or first_opt == "-23"):
        return "stateless"
    else:
        return "none"

def comm_input_output(options, stdin, stdout):
    first_opt = format_arg_chars(options[0])
    if(first_opt == "-13"):
        input_opt = format_arg_chars(options[2])
        if(input_opt == "-"):
            in_stream = ["stdin"]
            opt_indices = [("option", i) for i in range(len(options))]
        else:
            in_stream = [("option", 2)]
            opt_indices = [("option", 0), ("option", 1)]
        return (in_stream, ["stdout"], opt_indices)
    elif (first_opt == "-23"):
        input_opt = format_arg_chars(options[1])
        if(input_opt == "-"):
            in_stream = ["stdin"]
            opt_indices = [("option", i) for i in range(len(options))]
        else:
            in_stream = [("option", 1)]
            opt_indices = [("option", 0), ("option", 2)]
        return (in_stream, ["stdout"], opt_indices)
    else:
        assert(false)

## WARNING: This is not complete!! It doesn't handle - for stdin, or
## directories, etc...
##
## TODO: Make this complete
def diff_input_output(options, stdin, stdout):
    in_stream = [("option", i) for i, option in enumerate(options)
                 if not format_arg_chars(option).startswith('-')]
    opt_indices = [("option", i) for i, option in enumerate(options)
                   if format_arg_chars(option).startswith('-')]
    return (in_stream, ["stdout"], opt_indices)

def make_split_files(in_fid, fan_out, batch_size, fileIdGen):
    assert(fan_out > 1)
    ## Generate the split file ids
    out_fids = [fileIdGen.next_file_id() for i in range(fan_out)]
    split_commands = []
    curr = in_fid
    out_i = 0
    while (out_i + 2 < len(out_fids)):
        temp_fid = fileIdGen.next_file_id()
        split_com = make_split_file(curr, [out_fids[out_i], temp_fid], batch_size)
        split_commands.append(split_com)

        curr = temp_fid
        out_i += 1

    ## The final 2 children of out_fid
    split_com = make_split_file(curr, out_fids[out_i:(out_i+2)], batch_size)
    split_commands.append(split_com)
    return split_commands, out_fids

## TODO: Make a proper splitter subclass of Node
def make_split_file(in_fid, out_fids, batch_size):
    assert(len(out_fids) == 2)
    ## TODO: Call split_file recursively when we want to split a file
    ## more than two times

    ## TODO: I probably have to give the file names as options to the command to.
    options = [string_to_argument(str(batch_size))] + out_fids
    opt_indices = [("option", 0)]
    command = Command(None, # TODO: Make a proper AST
                      string_to_argument("split_file"),
                      options,
                      ["stdin"],
                      [("option", 1), ("option", 2)],
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

    def serialize_as_JSON(self):
        output_json = {}
        nodes = {}
        all_file_ids = []
        for i, node in enumerate(self.nodes):

            ## Gather all pipe names so that they are generated in the
            ## backend.
            input_pipes = [fid.serialize()
                           for fid in node.get_input_file_ids()
                           if fid.resource is None]
            output_pipes = [fid.serialize()
                            for fid in node.get_output_file_ids()
                            if fid.resource is None]
            all_file_ids += input_pipes + output_pipes

            ## Find the stdin and stdout files of nodes so that the
            ## backend can make the necessary redirections.
            if ("stdin" in node.in_stream):
                stdin_input_pipes = [Find(node.stdin)]
            else:
                stdin_input_pipes = []

            if ("stdout" in node.out_stream):
                stdout_output_pipes = [Find(node.stdout)]
            else:
                stdout_output_pipes = []

            ## TODO: Check if leaving the stdin and stdout pipes could
            ## lead to any problem.
            node_json = {}
            node_json["in"] = stdin_input_pipes
            node_json["out"] = stdout_output_pipes
            node_json["command"] = node.serialize()
            nodes[str(i)] = node_json

        all_file_ids = list(set(all_file_ids))
        output_json["fids"] = all_file_ids
        output_json["nodes"] = nodes
        flat_stdin = [Find(fid) for fid in self.stdin]
        output_json["in"] = [fid.serialize() for fid in flat_stdin]
        flat_stdout = [Find(fid) for fid in self.stdout]
        output_json["out"] = [fid.serialize() for fid in flat_stdout]
        return output_json

    def serialize_as_JSON_string(self):
        output_json = self.serialize_as_JSON()
        return json.dumps(output_json, sort_keys=True, indent=4)

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
            for other_node in self.nodes:
                ## Note: What if other_node == node?
                if (not incoming_fid.find_fid_list(other_node.get_output_file_ids()) is None):
                    previous_nodes_and_edges.append((other_node, incoming_fid))
        return previous_nodes_and_edges

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

