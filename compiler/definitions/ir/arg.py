from ir_utils import *

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
    
    def to_ast(self):
        return self.arg_char_list

