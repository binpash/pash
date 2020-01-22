import copy
import json
from union_find import *
from util import *

### TODO: Move this somewhere else
stateless_commands = ["cat", "tr", "grep", "col",
                      "groff", # not clear
                      "sed", # not always
                      "cut",
                      "gunzip", # file stateless
                      "xargs"] # I am not sure if all xargs are stateless

pure_commands = ["sort", "wc", "uniq", "bigrams_aux", "alt_bigrams_aux"]
parallelizable_pure_commands = ["sort", "bigrams_aux", "alt_bigrams_aux"]


### Utils

## This function gets a key and a value from a dictionary that only
## has one key
def get_kv(dic):
    dic_items = list(dic.items())
    assert(len(dic_items) == 1)
    return dic_items[0]

def format_arg_chars(arg_chars):
    chars = [format_arg_char(arg_char) for arg_char in arg_chars]
    return "".join(chars)

def format_arg_char(arg_char):
    key, val = get_kv(arg_char)
    if (key == 'C'):
        return str(chr(val))
    elif (key == 'B'):
        # The $() is just for illustration. This is backticks
        return '$({})'.format(val)
    elif (key == 'Q'):
        return '"{}"'.format(format_arg_chars(val))
    elif (key == 'V'):
        return '${{{}}}'.format(val[2])
    elif (key == 'E'):
        return '{}'.format(chr(val))
    else:
        ## TODO: Make this correct
        return key

def string_to_argument(string):
    return [char_to_arg_char(char) for char in string]

def make_argument(option):
    if(isinstance(option, FileId)):
        ## This is how commands are initialized when duplicating
        return option
    else:
        ## This is how commands are initialized in the AST
        return Arg(option)


## FIXME: This is certainly not complete. It is used to generate the
## AST for the call to the distributed planner. It only handles simple
## characters
def char_to_arg_char(char):
    return { 'C' : ord(char) }

## TODO: Resources should probably be more elaborate than just a
## string and a line range. They could be URLs, and possibly other things.
class Resource:
    def __init__(self, uri, range = [0, "inf"]):
        self.uri = uri
        self.range = range

    def __repr__(self):
        if self.range[0] == 0 and self.range[1] == "inf":
            output = str(self.uri)
        else:
            output = "{}[{}:{}]".format(self.uri, self.range[0], self.range[1])
        return output

    ## This function splits a "splittable" resource (at the moment
    ## only files are resources so they are all splittable)
    def split_resource(self, times, batch_size):
        assert(times > 0 and batch_size > 0)
        init = self.range[0]
        resources = [Resource(self.uri, [i * batch_size + init, (i+1) * batch_size + init])
                     for i in range(times)]
        ## Set the end to be equal to the old end
        resources[-1].range[1] = self.range[1]
        return resources

    ## This returns the length of the resource or None if it is
    ## infinite.
    def get_length(self):
        if(self.range[1] == "inf"):
            return None
        else:
            return (self.range[1] - self.range[0])

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


## Note: The NULL ident is considered to be the default unknown file id
##
## TODO: WARNING: We have to make sure that a resource in our IR can
## be uniquely determined given a relative or absolute path. Actually,
## we need to make sure that expanding any variable/string in our IR,
## will always return the same result.
##
## WARNING: At the moment it is not clear what resources are saved in
## the Find(self) and in self. This might create problems down the
## road.
##
## TODO: When doing union, I have to really make both file ids point
## to the same file.
class FileId:
    def __init__(self, ident, resource=None, children = [], max_length = None):
        self.ident = ident
        ## Initialize the parent
        MakeSet(self)
        self.resource=resource
        self.children = children
        ## Max length shows what is the maximum possible length that
        ## this file shows to. Its use is mostly to split intermediate
        ## streams.
        self.max_length = max_length

    def __repr__(self):
        ## Note: Outputs the parent of the union and not the file id
        ##       itself.
        if (self.resource is None):
            if(self.max_length is None):
                output = "#file{}".format(Find(self).ident)
            else:
                output = "#file{}[max:{}]".format(Find(self).ident, self.max_length)
        else:
            output = "#file{}({})".format(Find(self).ident, self.resource.__repr__())
        return output

    def serialize(self):
        if (self.resource is None):
            if(self.max_length is None):
                output = "#file{}".format(Find(self).ident)
            else:
                output = "#file{}[max:{}]".format(Find(self).ident, self.max_length)
        else:
            output = "{}".format(self.resource)
        return output

    ## Serialize as an option for the JSON serialization when sent to
    ## the backend. This is needed as options can either be files or
    ## arguments, and in each case there needs to be a different
    ## serialization procedure.
    def opt_serialize(self):
        return '"{}"'.format(self.serialize())

    ## TODO: Maybe this can be merged with serialize
    def pipe_name(self):
        assert(self.resource is None)
        output = '"#file{}"'.format(Find(self).ident)
        return output

    def set_resource(self, resource):
        ## The resource cannot be reset. A pointer can never point to
        ## more than one resource.
        assert(self.resource is None)
        self.resource = resource

    def get_resource(self):
        return self.resource

    ## TODO: We might need to reconstruct the parents from children,
    ## so we might have to add a parent field in file ids.
    def set_children(self, children):
        assert(self.children == [])
        self.children = children

    def get_children(self):
        return self.children

    def split_resource(self, times, batch_size, fileIdGen):
        if(not self.resource is None):
            ## This works as expected if the file points to a resource
            resources = self.resource.split_resource(times, batch_size)
            split_file_ids = [create_file_id_for_resource(resource, fileIdGen)
                              for resource in resources]
            self.set_children(split_file_ids)
        else:
            ## If the file doesn't point to a resource (meaning that
            ## it is an intermediate file), then we just restrict its
            ## max length, and the implementation should know to only
            ## pass so many lines in this file.
            split_file_ids = [create_split_file_id(batch_size, fileIdGen)
                              for i in range(times)]
            split_file_ids[-1].set_max_length(None)
            self.set_children(split_file_ids)

    ## This must be used by the implementation to only transfer
    ## max_length lines in this file.
    def set_max_length(self, max_length):
        self.max_length = max_length

    def toFileName(self, prefix):
        output = "{}_file{}".format(prefix, Find(self).ident)
        return output

    def isNull(self):
        return self.ident == "NULL"

    ## TODO: This union-find structure is very brittle. It could break
    ## very easily and it is difficult to reason about the
    ## files. Replace this with some other structure. Maybe a map that
    ## keeps the unifications of file identifiers. Make sure that when
    ## uniting two files, their children and resources are modified
    ## accordingly.
    def union(self, other):
        Union(self, other)
        my_resource = self.get_resource()
        other_resource = Find(other).get_resource()
        ## It shouldn't be the case that both resources are not NULL
        assert(my_resource is None or
               other_resource is None or
               my_resource == other_resource)

        if (my_resource is None):
            self.set_resource(other_resource)
        elif (other_resource is None):
            Find(other).set_resource(my_resource)

    def find_fid_list(self, fids):
        parent_fids = [Find(other_fid) for other_fid in fids]
        try:
            return parent_fids.index(Find(self))
        except ValueError:
            return None

    def get_ident(self):
        return self.ident

    def flatten(self):
        if(len(self.get_children()) > 0):
            return flatten_list([child.flatten() for child in self.get_children()])
        else:
            return [self]


class FileIdGen:
    def __init__(self, next = 0):
        self.next = next + 1

    def next_file_id(self):
        fileId = FileId(self.next)
        self.next += 1
        return fileId

## Question: Is this information adequate?
##
## TODO: What other information should a node of the IR contain?
## (other redirections possibly?).
##
## (LATER) TODO: Replace all the file references in IR nodes with new
## Identifiers that we make. IN order to do this, one has to be able
## to find these file arguments (based on the analysis that we will
## do).
##
## A node represents an abstract program that our system can
## distribute. At the moment, that is a program with one input and one
## output stream. Input and output streams are shown as a list of
## either options or standard channels (such as stdin, stdout,
## stderr).
##
## Nodes also have a category, which shows whether they can be
## parallelized on their input stream or not.
class Node:
    def __init__(self, ast, in_stream=[], out_stream=[],
                 category="none", stdin=None, stdout=None):
        self.ast = ast
        self.in_stream = in_stream
        self.out_stream = out_stream
        self.stdin = stdin
        self.stdout = stdout
        self.category = category

    def __repr__(self):
        output = "Node: \"{}\" in:{} out:{}".format(
            self.ast, self.stdin, self.stdout)
        return output

    ## These two commands return the flattened fileId list. Meaning
    ## that they return the children, if they exist.
    def get_flat_input_file_ids(self):
        return flatten_list([Find(file_id).flatten() for file_id in self.get_input_file_ids()])

    def get_flat_output_file_ids(self):
        return flatten_list([Find(file_id).flatten() for file_id in self.get_output_file_ids()])

    def get_input_file_ids(self):
        return [self.get_file_id(input_chunk) for input_chunk in self.in_stream]

    def get_output_file_ids(self):
        return [self.get_file_id(output_chunk) for output_chunk in self.out_stream]

    ## TODO: Rename
    def get_file_id(self, chunk):
        if (chunk == "stdout"):
            return self.stdout
        elif (chunk == "stdin"):
            return self.stdin
        elif (isinstance(chunk, tuple)
              and len(chunk) == 2
              and chunk[0] == "option"):
            ## If an option is asked, this node must be a command.
            assert(isinstance(self, Command))
            return self.options[chunk[1]]
        else:
            ## TODO: Complete this
            print(chunk)
            assert(False)

    ## TODO: Is there a way to abstract the behaviour of these two functions?
    def set_file_id(self, chunk, value):
        if (chunk == "stdout"):
            self.stdout = value
        elif (chunk == "stdin"):
            self.stdin = value
        elif (isinstance(chunk, tuple)
              and len(chunk) == 2
              and chunk[0] == "option"):
            ## If an option is asked, this node must be a command.
            assert(isinstance(self, Command))
            self.options[chunk[1]] = value
        else:
            ## TODO: Complete this
            print(chunk, value)
            assert(False)

    def find_file_id_in_in_stream(self, fileId):
        return self.find_file_id_in_stream(fileId, self.in_stream)

    def find_file_id_in_out_stream(self, fileId):
        return self.find_file_id_in_stream(fileId, self.out_stream)

    def find_file_id_in_stream(self, file_id, stream):
        index = 0
        for chunk in stream:
            chunk_file_id = Find(self.get_file_id(chunk))
            flat_file_ids = chunk_file_id.flatten()
            if(file_id in flat_file_ids):
                return index
            index += 1
        return None


## Commands are specific Nodes that can be parallelized if they are
## classified as stateless, etc...
class Command(Node):
    def __init__(self, ast, command, options, in_stream, out_stream,
                 opt_indices, category, stdin=None, stdout=None):
        super().__init__(ast, in_stream, out_stream, category, stdin, stdout)
        self.command = Arg(command)
        self.options = [make_argument(opt) for opt in options]
        self.opt_indices = opt_indices

    def __repr__(self):
        prefix = "Command"
        if (self.category == "stateless"):
            prefix = "Stateless"
        elif (self.category == "pure"):
            prefix = "Pure"
        # output = "{}: \"{}\" in:{} out:{} opts:{}".format(
        #     prefix, self.command, self.stdin, self.stdout, self.options)
        output = "{}: \"{}\" in:{} out:{}".format(
            prefix, self.command, self.get_flat_input_file_ids(),
            self.get_flat_output_file_ids())
        return output

    def serialize(self):
        all_opt_indices = [o_i[1] for o_i in (self.opt_indices + self.in_stream + self.out_stream)
                           if isinstance(o_i, tuple)]
        all_opt_indices.sort()
        options_string = " ".join([self.options[opt_i].opt_serialize() for opt_i in all_opt_indices])
        output = "{} {}".format(self.command, options_string)
        return output

    def get_non_file_options(self):
        return [self.options[i] for _, i in self.opt_indices]

    ## Get the file names of the outputs of the map commands. This
    ## differs if the command is stateless, pure that can be
    ## written as a map and a reduce, and a pure that can be
    ## written as a generalized map and reduce.
    def get_map_output_files(self, input_file_ids, fileIdGen):
        assert(self.category == "stateless" or self.is_pure_parallelizable())
        if(self.category == "stateless"):
            return [fileIdGen.next_file_id() for in_fid in input_file_ids]
        elif(self.is_pure_parallelizable()):
            return self.pure_get_map_output_files(input_file_ids, fileIdGen)
        else:
            print("Unreachable code reached :(")
            assert(False)
            ## This should be unreachable

    def pure_get_map_output_files(self, input_file_ids, fileIdGen):
        assert(self.is_pure_parallelizable())
        if(str(self.command) == "sort"):
            new_output_file_ids = [[fileIdGen.next_file_id()] for in_fid in input_file_ids]
        elif(str(self.command) == "bigrams_aux"):
            new_output_file_ids = [[fileIdGen.next_file_id() for i in range(BigramGMap.num_outputs)]
                                   for in_fid in input_file_ids]
        elif(str(self.command) == "alt_bigrams_aux"):
            new_output_file_ids = [[fileIdGen.next_file_id()] for in_fid in input_file_ids]
        else:
            print("Unreachable code reached :(")
            assert(False)
            ## This should be unreachable
        return new_output_file_ids

    def duplicate(self, new_output_file_ids, fileIdGen):
        assert(self.category == "stateless" or self.is_pure_parallelizable())
        if(self.category == "stateless"):
            return self.stateless_duplicate(new_output_file_ids)
        elif(self.is_pure_parallelizable()):
            return self.pure_duplicate(new_output_file_ids, fileIdGen)
        else:
            print("Unreachable code reached :(")
            assert(False)
            ## This should be unreachable

    def stateless_duplicate(self, output_file_ids):
        assert(self.category == "stateless")

        ## Attach the new output files as children of the node's
        ## output, because make_duplicate command requires that. Also,
        ## by doing that, the input of the next command now also has
        ## children (due to unification).
        out_edge_file_ids = self.get_output_file_ids()
        assert(len(out_edge_file_ids) == 1)
        out_edge_file_id = out_edge_file_ids[0]
        out_edge_file_id.set_children(output_file_ids)

        input_file_ids = self.get_flat_input_file_ids()

        in_out_file_ids = zip(input_file_ids, output_file_ids)

        new_commands = [self.make_duplicate_command(in_fid, out_fid) for in_fid, out_fid in in_out_file_ids]

        return new_commands

    def is_pure_parallelizable(self):

        ## TODO: Read from some file that contains information about
        ## commands instead of hardcoding
        return (self.category == "pure" and str(self.command) in parallelizable_pure_commands)

    def pure_duplicate(self, output_file_ids, fileIdGen):
        assert(self.is_pure_parallelizable())
        input_file_ids = self.get_flat_input_file_ids()

        in_out_file_ids = zip(input_file_ids, output_file_ids)

        simple_map_pure_commands = ["sort",
                                    "alt_bigrams_aux"]
        ## This is the category of all commands that don't need a
        ## special generalized map
        if(str(self.command) in simple_map_pure_commands):

            ## make_duplicate_command duplicates a node based on its
            ## output file ids, so we first need to assign them.
            intermediate_output_file_id = fileIdGen.next_file_id()
            new_output_file_ids = [fids[0] for fids in output_file_ids]
            intermediate_output_file_id.set_children(new_output_file_ids)
            self.stdout = intermediate_output_file_id

            new_commands = [self.make_duplicate_command(in_fid, out_fids[0])
                            for in_fid, out_fids in in_out_file_ids]
        elif(str(self.command) == "bigrams_aux"):
            new_commands = [BigramGMap([in_fid] + out_fids)
                            for in_fid, out_fids in in_out_file_ids]
        else:
            print("Unreachable code reached :(")
            assert(False)
            ## This should be unreachable
        return new_commands


    def make_duplicate_command(self, in_fid, out_fid):

        ## First find what does the new file identifier refer to
        ## (stdin, or some argument)
        new_in_stream_index = self.find_file_id_in_in_stream(in_fid)
        new_out_stream_index = self.find_file_id_in_out_stream(out_fid)
        new_in_stream = [self.in_stream[new_in_stream_index]]
        new_out_stream = [self.out_stream[new_out_stream_index]]

        ## TODO: Simplify the code below
        in_chunk = new_in_stream[0]
        if(in_chunk == "stdin"):
            new_stdin = in_fid
        else:
            ## Question: Is that valid?
            new_stdin = self.stdin

        if(new_out_stream[0] == "stdout"):
            new_stdout = out_fid
        else:
            ## Question: Is that valid?
            new_stdout = self.stdout

        new_options = self.options.copy()
        if(isinstance(in_chunk, tuple)
           and len(in_chunk) == 2
           and in_chunk[0] == "option"):
            new_options[in_chunk[1]] = in_fid

        ## TODO: I probably have to do the same with output options

        new_command = Command(None, # The ast is None
                              self.command,
                              new_options,
                              new_in_stream,
                              new_out_stream,
                              self.opt_indices,
                              self.category,
                              new_stdin,
                              new_stdout)
        ## Question: Is it valid setting stdin and stdout to the stdin
        ## and stdout of the current command?
        return new_command

###
### This part of the file contains special nodes. Either ones that
### will be later parsed by the DSL, or core nodes like cat, tee, etc.
###

## (Maybe) Extend this to also take flags as part of its input.
class Cat(Command):
    def __init__(self, file_ids):
        command = string_to_argument("cat")
        options = file_ids
        in_stream = [("option", i)  for i in range(len(file_ids))]
        out_stream = ["stdout"]
        opt_indices = []
        category = "stateless"
        ## TODO: Fill the AST
        ast = None
        super().__init__(ast, command, options, in_stream, out_stream,
                         opt_indices, category)

## These are the generalized map and reduce nodes for the
## `tee s2 | tail +2 | paste s2 -` function. At some point,
## these should be parsed from some configuration file, but for now
## we have them hardcoded for commands that we are interested in.
class BigramGMap(Command):
    num_outputs = 3
    def __init__(self, file_ids):
        command = string_to_argument("bigram_aux_map")
        options = [string_to_argument(fid.pipe_name()) for fid in file_ids]
        ## TODO: Generalize the model to arbitrarily many outputs to
        ## get rid of this hack, where inputs, outputs are just
        ## written as options.
        in_stream = []
        out_stream = []
        ## TODO: When we generalize the model to have arbitrary
        ## outputs we can have them also be outputs.
        opt_indices = [("option", 0), ("option", 1), ("option", 2), ("option", 3)]
        category = "stateless"
        super().__init__(None, command, options, in_stream, out_stream,
                         opt_indices, category)

class BigramGReduce(Command):
    def __init__(self, file_ids):
        command = string_to_argument("bigram_aux_reduce")
        options = self.make_options(file_ids)
        ## TODO: Generalize the model to arbitrarily many outputs to
        ## get rid of this hack, where inputs, outputs are just
        ## written as options.
        in_stream = []
        ## WARNING: This cannot be changes to be empty, because then
        ## it would not be correctly unified in the backend with the
        ## input of its downstream node.
        out_stream = [("option", 6)]
        ## TODO: When we generalize the model to have arbitrary
        ## outputs we can have them also be outputs.
        opt_indices = [("option", i) for i in range(0,9) if not i == 6]
        category = "pure"
        super().__init__(None, command, options, in_stream, out_stream,
                         opt_indices, category)

    ## TODO: Make this atrocity prettier
    def make_options(self, file_ids):
        options = []
        for i in range(0, 9):
            fid = file_ids[i]
            if(not i == 6):
                opt = string_to_argument(fid.pipe_name())
            else:
                opt = fid
            options.append(opt)
        return options

class AltBigramGReduce(Command):
    def __init__(self, file_ids):
        command = string_to_argument("alt_bigram_aux_reduce")
        options = self.make_options(file_ids)
        ## TODO: Generalize the model to arbitrarily many outputs to
        ## get rid of this hack, where inputs, outputs are just
        ## written as options.
        in_stream = []
        ## WARNING: This cannot be changes to be empty, because then
        ## it would not be correctly unified in the backend with the
        ## input of its downstream node.
        out_stream = [("option", 2)]
        ## TODO: When we generalize the model to have arbitrary
        ## outputs we can have them also be outputs.
        opt_indices = [("option", 0), ("option", 1)]
        category = "pure"
        super().__init__(None, command, options, in_stream, out_stream,
                         opt_indices, category)

    ## TODO: Make this atrocity prettier
    def make_options(self, file_ids):
        options = []
        for i in range(0, 3):
            fid = file_ids[i]
            if(not i == 2):
                opt = string_to_argument(fid.pipe_name())
            else:
                opt = fid
            options.append(opt)
        return options


class SortGReduce(Command):
    def __init__(self, file_ids):
        command = string_to_argument("sort")
        input_file_ids = file_ids[:-1]
        output_file_id = file_ids[-1]
        input_file_id_opts = [string_to_argument(fid.pipe_name()) for fid in input_file_ids]
        options = [string_to_argument("-m")] + input_file_id_opts
        ## Question: Can putting an empty input stream be a problem?
        ## If we do that, then this node will be a source in the
        ## dataflow graph. Ideally we want to generalize to arbitrary
        ## inputs. This also applies to the one above.
        in_stream = []
        out_stream = ["stdout"]
        ## At the moment none of the $IN1, $IN2 are considered inputs
        ##
        ## TODO: When we generalize the model to have arbitrary
        ## inputs we can have them also be inputs.
        opt_indices = [("option", 0), ("option", 1), ("option", 2)]
        category = "pure"
        super().__init__(None, command, options, in_stream, out_stream,
                         opt_indices, category, stdout=output_file_id)


def create_command_assign_file_identifiers(ast, fileIdGen, command, options, stdin=None, stdout=None):
    in_stream, out_stream, opt_indices = find_command_input_output(command, options, stdin, stdout)
    category = find_command_category(command, options)
    command = Command(ast, command, options, in_stream, out_stream,
                      opt_indices, category, stdin, stdout)

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


class Arg:
    def __init__(self, arg_char_list):
        if(isinstance(arg_char_list, Arg)):
           ## TODO: We might need to copy here?
           self.arg_char_list = arg_char_list.arg_char_list
        else:
           self.arg_char_list = arg_char_list

    def __repr__(self):
        return format_arg_chars(self.arg_char_list)

    def opt_serialize(self):
        return self.__repr__()

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


    if (command_string in stateless_commands):
        return "stateless"
    elif (command_string in pure_commands):
        return "pure"
    elif (command_string == "comm"):
        return is_comm_pure(options)
    else:
        return "none"

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

def make_split_files(in_fid, out_fid, batch_size, fileIdGen):
    assert(len(out_fid.children) >= 2)
    split_commands = []
    curr = in_fid
    out_i = 0
    while (out_i + 2 < len(out_fid.children)):
        temp_fid = fileIdGen.next_file_id()
        temp_out_fid = fileIdGen.next_file_id()
        temp_out_fid.set_children([out_fid.children[out_i], temp_fid])
        split_com = make_split_file(curr, temp_out_fid, batch_size)
        split_commands.append(split_com)

        curr = temp_fid
        out_i += 1

    ## The final 2 children of out_fid
    final_out_fid = fileIdGen.next_file_id()
    final_out_fid.set_children(out_fid.children[out_i:(out_i+2)])
    split_com = make_split_file(curr, final_out_fid, batch_size)
    split_commands.append(split_com)
    return split_commands

## TODO: Make a proper splitter subclass of Node
def make_split_file(in_fid, out_fid, batch_size):
    assert(len(out_fid.children) == 2)
    ## TODO: Call split_file recursively when we want to split a file
    ## more than two times

    ## TODO: I probably have to give the file names as options to the command to.
    options = [string_to_argument(str(batch_size))]
    opt_indices = [("option", i) for i in range(len(options))]
    command = Command(None, # TODO: Make a proper AST
                      string_to_argument("split_file"),
                      options,
                      ["stdin"],
                      ["stdout"],
                      opt_indices,
                      None, # TODO: Category?
                      in_fid,
                      out_fid)
    return command

## This function gets a file identifier and returns the maximum among
## its, and its parents identifier (parent regarding Union Find)
def get_larger_file_id_ident(file_id):
    my_ident = file_id.get_ident()
    find_ident = Find(file_id).get_ident()
    return max(my_ident, find_ident)

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
    def __init__(self, nodes, stdin = None, stdout = None):
        self.nodes = nodes
        self.edges = {}
        if(stdin is None):
            self.stdin = FileId("NULL")
        else:
            self.stdin = stdin
        if(stdout is None):
            self.stdout = FileId("NULL")
        else:
            self.stdout = stdout

    def __repr__(self):
        output = "(|-{} IR: {} {}-|)".format(self.stdin.flatten(), self.nodes, self.stdout.flatten())
        return output

    def serialize(self):
        output = "Nodes:\n"
        all_file_ids = ""
        for i, node in enumerate(self.nodes):
            serialized_input_file_ids = " ".join([fid.serialize()
                                                  for fid in node.get_flat_input_file_ids()])
            serialized_output_file_ids = " ".join([fid.serialize()
                                                   for fid in node.get_flat_output_file_ids()])
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
                           for fid in node.get_flat_input_file_ids()
                           if fid.resource is None]
            output_pipes = [fid.serialize()
                            for fid in node.get_flat_output_file_ids()
                            if fid.resource is None]
            all_file_ids += input_pipes + output_pipes

            ## Find the stdin and stdout files of nodes so that the
            ## backend can make the necessary redirections.
            if ("stdin" in node.in_stream):
                stdin_input_pipes = Find(node.stdin).flatten()
            else:
                stdin_input_pipes = []

            if ("stdout" in node.out_stream):
                stdout_output_pipes = Find(node.stdout).flatten()
            else:
                stdout_output_pipes = []

            node_json = {}
            node_json["in"] = stdin_input_pipes
            node_json["out"] = stdout_output_pipes
            node_json["command"] = node.serialize()
            nodes[str(i)] = node_json

        all_file_ids = list(set(all_file_ids))
        output_json["fids"] = all_file_ids
        output_json["nodes"] = nodes
        flat_stdin = flatten_list([Find(self.stdin).flatten()])
        output_json["in"] = [fid.serialize() for fid in flat_stdin]
        flat_stdout = flatten_list([Find(self.stdout).flatten()])
        output_json["out"] = [fid.serialize() for fid in flat_stdout]
        return output_json

    def serialize_as_JSON_string(self):
        output_json = self.serialize_as_JSON()
        return json.dumps(output_json, sort_keys=True, indent=4)

    def set_ast(self, ast):
        self.ast = ast

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
        assert(not self.stdout.isNull())
        assert(not other.stdin.isNull())
        self.stdout.union(other.stdin)
        self.stdout = other.stdout

        ## Note: The ast is not extensible, and thus should be
        ## invalidated if an operation happens on the IR
        self.ast = None

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

    def get_next_nodes(self, node):
        next_nodes = []
        for outgoing_fid in node.get_output_file_ids():
            for other_node in self.nodes:
                ## Note: What if other_node == node?
                if (not outgoing_fid.find_fid_list(other_node.get_input_file_ids()) is None):
                    next_nodes.append(other_node)
        return next_nodes

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
                max_id = max(get_larger_file_id_ident(file_id), max_id)
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
                not self.stdin.isNull() and
                not self.stdout.isNull())

