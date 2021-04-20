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

    def raw_str(self):
        return format_arg_to_raw_chars(self.arg_char_list)