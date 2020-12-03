from json import JSONEncoder

from definitions.ast_node_c import *
from definitions.no_match_exception import *
from ir_utils import *


## TODO: Create subclasses for all different types of AstNodes
class AstNode:
    # create an AstNode object from an ast object as parsed by libdash
    def __init__(self, ast_object):
        try:
            self.construct = AstNodeConstructor(ast_object[0])
            self.parse_args(ast_object[1])
        except ValueError as no_matching_construct:
            raise NoMatchException('{} is not a construct we can handle'.format(ast_object[0]))

    def parse_args(self, args):
        if self.construct is AstNodeConstructor.PIPE:
            self.is_background = args[0]
            self.items = args[1]
        elif self.construct is AstNodeConstructor.COMMAND:
            self.line_number = args[0]
            self.assignments = args[1]
            self.arguments = args[2]
            self.redir_list = args[3]
            return
        elif self.construct in [AstNodeConstructor.AND, AstNodeConstructor.OR, AstNodeConstructor.SEMI]:
            self.left_operand = args[0]
            self.right_operand = args[1]
        elif self.construct in [AstNodeConstructor.REDIR, AstNodeConstructor.SUBSHELL, AstNodeConstructor.BACKGROUND]:
            self.line_number = args[0]
            # TODO maybe pick a better name?
            self.node = args[1]
            self.redir_list = args[2]
        elif self.construct is AstNodeConstructor.DEFUN:
            self.line_number = args[0]
            self.name = args[1]
            self.body = args[2]
        elif self.construct is AstNodeConstructor.FOR:
            self.line_number = args[0]
            self.argument = args[1]
            self.body = args[2]
            self.variable = args[3]
        else:
            raise ValueError()
    
    def __repr__(self):
        if self.construct is AstNodeConstructor.PIPE:
            if (self.is_background):
                return "Background Pipe: {}".format(self.items)    
            else:
                return "Pipe: {}".format(self.items)
        elif self.construct is AstNodeConstructor.COMMAND:
            output = "Command: {}".format(self.arguments)
            if(len(self.assignments) > 0):
                output += ", ass[{}]".format(self.assignments)
            if(len(self.redir_list) > 0):
                output += ", reds[{}]".format(self.redir_list)
            return output
        elif self.construct is AstNodeConstructor.FOR:
            output = "for {} in {}; do ({})".format(self.variable, self.argument, self.body)
            return output
        return NotImplemented 

    def check(self, **kwargs):
        # user-supplied custom checks
        for key, value in kwargs.items():
            try:
                assert(value())
            except Exception as exc:
                print("check for {} construct failed at key {}".format(self.construct, key))
                raise exc
    
    def json_serialize(self):
        if self.construct is AstNodeConstructor.FOR:
            json_output = make_kv(self.construct.value,
                           [self.line_number,
                            self.argument,
                            self.body,
                            self.variable])
        elif self.construct is AstNodeConstructor.COMMAND:
            json_output = make_kv(self.construct.value,
                                  [self.line_number,
                                   self.assignments,
                                   self.arguments,
                                   self.redir_list])
        elif self.construct is AstNodeConstructor.PIPE:
            json_output = make_kv(self.construct.value,
                                  [self.is_background,
                                   self.items])
        else:
            print(self)
            json_output = NotImplemented
        return json_output

class CustomJSONEncoder(JSONEncoder):
    def default(self, obj):
        if isinstance(obj, AstNode):
            return obj.json_serialize()
        # Let the base class default method raise the TypeError
        return JSONEncoder.default(self, obj)