import copy
from union_find import *
### Utils

## TODO: Move to another file
def flatten(l):
    return [item for sublist in l for item in sublist]

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
        return '${}'.format(val[2])
    else:
        ## TODO: Make this correct
        return key

def string_to_argument(string):
    return [char_to_arg_char(char) for char in string]

## FIXME: This is certainly not complete. It is used to generate the
## AST for the call to the distributed planner. It only handles simple
## characters
def char_to_arg_char(char):
    return { 'C' : ord(char) }
    
## Note: The NULL ident is considered to be the default unknown file id
##
## TODO: WARNING: We have to make sure that a resource in our IR can
## be uniquely determined given a relative or absolute path. Actually,
## we need to make sure that expanding any variable/string in our IR,
## will always return the same result.
class FileId:
    def __init__(self, ident, resource=None):
        self.ident = ident
        ## Initialize the parent
        MakeSet(self)
        self.resource=resource

    def __repr__(self):
        ## Note: Outputs the parent of the union and not the file id
        ##       itself.
        if (self.resource is None):
            output = "#file{}".format(Find(self).ident)
        else:
            output = "#file{}({})".format(Find(self).ident, self.resource.__repr__())
        return output

    def set_resource(self, resource):
        ## The resource cannot be reset. A pointer can never point to
        ## more than one resource.
        assert(self.resource is None)
        self.resource = resource

    def get_resource(self):
        return self.resource
    
    def toFileName(self, prefix):
        output = "{}_file{}".format(prefix, Find(self).ident)
        return output
    
    def isNull(self):
        return self.ident == "NULL"

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
        
    
class FileIdGen:
    def __init__(self):
        self.next = 0

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
    

## Commands are specific Nodes that can be parallelized if they are
## classified as stateless, etc...
class Command(Node):
    def __init__(self, ast, command, options, in_stream, out_stream,
                 category, stdin=None, stdout=None):
        super().__init__(ast, in_stream, out_stream, category, stdin, stdout)
        self.command = Arg(command)
        self.options = [Arg(opt) for opt in options]
        
    def __repr__(self):
        prefix = "Command"
        if (self.category == "stateless"):
            prefix = "Stateless"
        output = "{}: \"{}\" in:{} out:{} opts:{}".format(
            prefix, self.command, self.stdin, self.stdout, self.options)
        return output
    
def create_command_assign_file_identifiers(ast, fileIdGen, command, options, stdin=None, stdout=None):
    in_stream, out_stream = find_command_input_output(command, options, stdin, stdout)
    category = find_command_category(command, options)
    command = Command(ast, command, options, in_stream, out_stream, category, stdin, stdout)

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
    fid = command.get_file_id(opt_or_channel)
    ## If the file is not a FileId, then it is some argument. We
    ## create a file identifier, and replace it with that, and
    ## make sure that the file identifier points to the argument.
    if (not isinstance(fid, FileId)):
        real_fid = fileIdGen.next_file_id()
        real_fid.set_resource(fid)
        return real_fid
    else:
        return fid


class Arg:
    def __init__(self, arg_char_list):
        self.arg_char_list = arg_char_list

    def __repr__(self):
        return format_arg_chars(self.arg_char_list)
    
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
        return (input_stream, ["stdout"])
    else:
        return (["stdin"], ["stdout"])


## This functions finds and returns a string representing the command category
def find_command_category(command, options):
    command_string = format_arg_chars(command)
    # print("Command to categorize:", command_string)

    assert(isinstance(command_string, str))

    ## TODO: Make a proper search that returns the command category
    print(" -- Warning: Category for: {} is hardcoded and possibly wrong".format(command_string))
    
    if (command_string == "cat"):
        return "stateless"
    else:
        return "none"

    
    
## Note: This might need more information. E.g. all the file
## descriptors of the IR, and in general any other local information
## that might be relevant.
class IR:
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
        output = "(|-{} IR: {} {}-|)".format(self.stdin, self.nodes, self.stdout)
        return output

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
                if (Find(incoming_fid) in
                    [Find(out_fid) for out_fid in other_node.get_output_file_ids()]):
                    return True
        return False
    
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

