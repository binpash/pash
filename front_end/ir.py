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

def format_arg_char(arg_char):
    key, val = get_kv(arg_char)
    if (key == 'C'):
        return str(chr(val))
    elif (key == 'B'):
        # The $() is just for illustration. This is backticks
        return '$({})'.format(val)
    elif (key == 'Q'):
        return '"{}"'.format("".join([format_arg_char(c) for c in val]))
    elif (key == 'V'):
        return '${}'.format(val[2])
    else:
        ## TODO: Make this correct
        return key

## Note: The NULL ident is considered to be the default unknown file id
class FileId:
    def __init__(self, ident):
        self.ident = ident
        ## Initialize the parent
        MakeSet(self)

    def __repr__(self):
        ## Note: Outputs the parent of the union and not the file id
        ##       itself.
        output = "#file{}".format(Find(self).ident)
        # output = "#file{}".format(self.ident)
        return output

    def isNull(self):
        return self.ident == "NULL"
    
class FileIdGen:
    def __init__(self):
        self.next = 0

    def next_file_id(self):
        fileId = FileId(self.next)
        self.next += 1
        return fileId
    
## Question: Is this information adequate?
##
## TODO: This is not really a command. It is a node in the IR. Rename
## this class to NodeIR. Its fields should be an AST, stdin, stdout
## (and other redirections possibly?).
##
## (LATER) TODO: Replace all the file references in IR nodes with new
## Identifiers that we make. IN order to do this, one has to be able
## to find these file arguments (based on the analysis that we will
## do).
class Command:
    def __init__(self, command, stdin=None, stdout=None, options=None):
        self.command = Arg(command)
        self.stdin = stdin
        self.stdout = stdout
        self.options = [Arg(opt) for opt in options]
        
    def __repr__(self):
        output = "Command: \"{}\" in:{} out:{} opts:{}".format(
            self.command, self.stdin, self.stdout, self.options)
        return output

class Arg:
    def __init__(self, arg_char_list):
        self.arg_char_list = arg_char_list

    def __repr__(self):
        chars = [format_arg_char(arg_char) for arg_char in self.arg_char_list]
        return "".join(chars)

    
## Note: This might need more information. E.g. all the file
## descriptors of the IR, and in general any other local information
## that might be relevant.
class IR:
    def __init__(self, nodes, stdin = None, stdout = None):
        self.nodes = nodes
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
        ## stdout of the second.
        ##
        ## Question: What happens if one of them is NULL. This
        ##           shouldn't be the case after we check that
        ##           both self and other are not empty.
        assert(not self.stdout.isNull())
        assert(not other.stdin.isNull())
        Union(self.stdout, other.stdin)
        self.stdout = other.stdout

        ## Note: The ast is not extensible, and thus should be
        ## invalidated if an operation happens on the IR
        self.ast = None
            
    ## Note: We assume that no nodes is an adequate condition to check
    ##       emptiness.
    def empty(self):
        return (len(self.nodes) == 0)

    ## This function checks whether an IR is valid -- that is, if it
    ## has at least one node, and stdin, stdout set to some non-null
    ## file identifiers.
    def valid(self):
        return (len(self.nodes) > 0 and
                not self.stdin.isNull() and
                not self.stdin.isNull())
