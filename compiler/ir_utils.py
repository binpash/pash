### Utils

def get_command_from_definition(command_definition):
    if 'command' in command_definition:
        return command_definition['command']

    print('Possible issue with definition file: Missing command in command definition {}'.format(command_definition))
    return ''


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
        ## TODO: This is not right. I think the main reason for the
        ## problems is the differences between bash and the posix
        ## standard.
        # print(" -- escape-debug -- ", val, chr(val))
        non_escape_chars = [92, # \
                            61, # =
                            91, # [
                            93, # ]
                            45, # -
                            58, # :
                            42] # *
        if(val in non_escape_chars):
            return '{}'.format(chr(val))
        else:
            return '\{}'.format(chr(val))
    else:
        ## TODO: Make this correct
        return key

## This function gets a key and a value from the ast json format
def get_kv(dic):
    return (dic[0], dic[1])

def make_kv(key, val):
    return [key, val]


def string_to_argument(string):
    return [char_to_arg_char(char) for char in string]

## FIXME: This is certainly not complete. It is used to generate the
## AST for the call to the distributed planner. It only handles simple
## characters
def char_to_arg_char(char):
    return ['C' , ord(char)]

