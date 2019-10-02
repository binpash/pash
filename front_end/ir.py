
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
    else:
        ## TODO: Make this correct
        return key

## Question: Is this information adequate?
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

## Note: This might need more information to simplify the
## implementation (like the IRs stdin, stdout)
class IR:
    def __init__(self, nodes = []):
        self.nodes = nodes

    def __repr__(self):
        output = "(| IR: {} |)".format(self.nodes)
        return output

    ## TODO: There have to be more complex methods to combine two IRs
    def append(self, other):
        self.nodes += other.nodes
        return self

    def empty(self):
        return (len(self.nodes) == 0)
