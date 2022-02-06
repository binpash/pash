from ir_utils import *
from util import *

class Arg:
    def __init__(self, arg_char_list):
        if(isinstance(arg_char_list, Arg)):
            ## TODO: Make sure that this does not happen using an assertion.
            ## TODO: We might need to copy here?
            self.arg_char_list = arg_char_list.arg_char_list
        else:
           self.arg_char_list = arg_char_list

    def __repr__(self):
        return format_arg_chars(self.arg_char_list)

    def __eq__(self, other):
        if(isinstance(other, Arg)):
            return self.arg_char_list == other.arg_char_list
        log("Warning: Comparing Arg:", self, "with a non Arg argument:", other, "of type:", type(other))
        return False

    def opt_serialize(self):
        return self.__repr__()
    
    def to_ast(self):
        return self.arg_char_list

    def concatenate(self, other):
        space = [['C', 32]]
        self.arg_char_list.extend(space)
        self.arg_char_list.extend(other.arg_char_list)

    def wrap_in_quotes(self):
        quote = Arg(string_to_argument("\"")).arg_char_list[0]
        escape = Arg(string_to_argument("\\")).arg_char_list[0]
        space = Arg(string_to_argument(" ")).arg_char_list[0]

        new_char_list = [quote]
        i = 0
        while i < len(self.arg_char_list):
            char = self.arg_char_list[i]
            if char == escape: 
                # This condition is useful in cases where we have repeated escape characters
                # example: "\\\""
                new_char_list.append(self.arg_char_list[i])
                i += 1 # #skip the next character
                new_char_list.append(self.arg_char_list[i])
            elif char == quote:
                new_char_list.extend([escape, quote])
            elif char[0] == 'Q':
                # Note: parsing and unparsing here because just removing the 'Q' key without reparsing
                # lead to wrong character escapes
                inside_quotes = string_to_argument(format_arg_chars([char]))[1:-1] # Should this be recursive?
                new_char_list.extend([escape, quote])
                new_char_list.extend(inside_quotes)
                new_char_list.extend([escape, quote])
            else:
                new_char_list.append(char)
            i += 1
        new_char_list.extend([quote]) # ending quote
        return Arg(new_char_list)
