from __future__ import annotations
from shell_ast.ast_util import *
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

    @staticmethod
    def string_to_arg(string) -> Arg:
        return Arg(string_to_argument(string))


