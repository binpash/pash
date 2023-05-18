from json import JSONEncoder

from shell_ast.ast_node_c import *
from definitions.no_match_exception import *
from util import *


class AstNode:
    NodeName = 'None'
    # create an AstNode object from an ast object as parsed by libdash
    # def __init__(self, ast_object):
    #     try:
    #         self.construct = AstNodeConstructor(ast_object[0])
    #         self.parse_args(ast_object[1])
    #     except ValueError as no_matching_construct:
    #         raise NoMatchException('{} is not a construct we can handle'.format(ast_object[0]))

    # def parse_args(self, args):
    #     if self.construct is AstNodeConstructor.PIPE:
    #         self.is_background = args[0]
    #         self.items = args[1]
    #     elif self.construct is AstNodeConstructor.COMMAND:
    #         self.line_number = args[0]
    #         self.assignments = args[1]
    #         self.arguments = args[2]
    #         self.redir_list = args[3]
    #     elif self.construct is AstNodeConstructor.SUBSHELL:
    #         self.line_number = args[0]
    #         self.body = args[1]
    #         self.redir_list = args[2]
    #     elif self.construct in [AstNodeConstructor.AND, AstNodeConstructor.OR, AstNodeConstructor.SEMI]:
    #         self.left_operand = args[0]
    #         self.right_operand = args[1]
    #     elif self.construct is AstNodeConstructor.NOT:
    #         self.body = args
    #     elif self.construct in [AstNodeConstructor.REDIR, AstNodeConstructor.BACKGROUND]:
    #         self.line_number = args[0]
    #         # TODO maybe pick a better name?
    #         self.node = args[1]
    #         self.redir_list = args[2]
    #     elif self.construct is AstNodeConstructor.DEFUN:
    #         self.line_number = args[0]
    #         self.name = args[1]
    #         self.body = args[2]
    #     elif self.construct is AstNodeConstructor.FOR:
    #         self.line_number = args[0]
    #         self.argument = args[1]
    #         self.body = args[2]
    #         self.variable = args[3]
    #     elif self.construct is AstNodeConstructor.WHILE:
    #         self.test = args[0]
    #         self.body = args[1]
    #     elif self.construct is AstNodeConstructor.IF:
    #         self.cond = args[0]
    #         self.then_b = args[1]
    #         self.else_b = args[2]
    #     elif self.construct is AstNodeConstructor.CASE:
    #         self.line_number = args[0]
    #         self.argument = args[1]
    #         self.cases = args[2]
    #     else:
    #         raise ValueError()
    
    # def __repr__(self):
    #     if self.construct is AstNodeConstructor.FOR:
    #         output = "for {} in {}; do ({})".format(self.variable, self.argument, self.body)
    #         return output
    #     log(self.construct)
    #     return NotImplemented 

    # def check(self, **kwargs):
    #     # user-supplied custom checks
    #     for key, value in kwargs.items():
    #         try:
    #             assert(value())
    #         except Exception as exc:
    #             log("check for {} construct failed at key {}".format(self.construct, key))
    #             raise exc
    
    # def json_serialize(self):
    #     if self.construct is AstNodeConstructor.FOR:
    #         json_output = make_kv(self.construct.value,
    #                        [self.line_number,
    #                         self.argument,
    #                         self.body,
    #                         self.variable])
    #     elif self.construct is AstNodeConstructor.WHILE:
    #         json_output = make_kv(self.construct.value,
    #                        [self.test,
    #                         self.body])
    #     # elif self.construct is AstNodeConstructor.COMMAND:
    #     #     json_output = make_kv(self.construct.value,
    #     #                           [self.line_number,
    #     #                            self.assignments,
    #     #                            self.arguments,
    #     #                            self.redir_list])
    #     # elif self.construct is AstNodeConstructor.REDIR:
    #     #     json_output = make_kv(self.construct.value,
    #     #                           [self.line_number,
    #     #                            self.node,
    #     #                            self.redir_list])
    #     # elif self.construct is AstNodeConstructor.BACKGROUND:
    #     #     json_output = make_kv(self.construct.value,
    #     #                           [self.line_number,
    #     #                            self.node,
    #     #                            self.redir_list])
    #     # elif self.construct is AstNodeConstructor.SUBSHELL:
    #     #     json_output = make_kv(self.construct.value,
    #     #                           [self.line_number,
    #     #                            self.body,
    #     #                            self.redir_list])
    #     # elif self.construct is AstNodeConstructor.PIPE:
    #     #     json_output = make_kv(self.construct.value,
    #     #                           [self.is_background,
    #     #                            self.items])
    #     # elif self.construct is AstNodeConstructor.DEFUN:
    #     #     json_output = make_kv(self.construct.value,
    #     #                           [self.line_number,
    #     #                            self.name,
    #     #                            self.body])
    #     # elif self.construct is AstNodeConstructor.IF:
    #     #     json_output = make_kv(self.construct.value,
    #     #                           [self.cond,
    #     #                            self.then_b,
    #     #                            self.else_b])
    #     # elif self.construct is AstNodeConstructor.SEMI:
    #     #     json_output = make_kv(self.construct.value,
    #     #                           [self.left_operand,
    #     #                            self.right_operand])
    #     # elif self.construct is AstNodeConstructor.OR:
    #     #     json_output = make_kv(self.construct.value,
    #     #                           [self.left_operand,
    #     #                            self.right_operand])
    #     # elif self.construct is AstNodeConstructor.AND:
    #     #     json_output = make_kv(self.construct.value,
    #     #                           [self.left_operand,
    #     #                            self.right_operand])
    #     # elif self.construct is AstNodeConstructor.NOT:
    #     #     json_output = make_kv(self.construct.value,
    #     #                           self.body)
    #     # elif self.construct is AstNodeConstructor.CASE:
    #     #     json_output = make_kv(self.construct.value,
    #     #                           [self.line_number,
    #     #                            self.argument,
    #     #                            self.cases])
    #     else:
    #         log("Not implemented serialization", self)
    #         json_output = NotImplemented
    #     return json_output

class CustomJSONEncoder(JSONEncoder):
    def default(self, obj):
        if isinstance(obj, AstNode):
            return obj.json_serialize()
        # Let the base class default method raise the TypeError
        return JSONEncoder.default(self, obj)

## TODO: Move serialize and repr to the subclasses below

## TODO: Determine all field types

class PipeNode(AstNode):
    NodeName = 'Pipe'
    is_background: bool
    items: object ## TODO: Find type

    def __init__(self, is_background, items):
        self.is_background = is_background
        self.items = items

    def __repr__(self):
        if (self.is_background):
            return "Background Pipe: {}".format(self.items)    
        else:
            return "Pipe: {}".format(self.items)
        
    def json_serialize(self):
        json_output = make_kv(PipeNode.NodeName,
                              [self.is_background,
                               self.items])
        return json_output

class CommandNode(AstNode):
    NodeName = 'Command'
    line_number: int
    assignments: object
    arguments: object
    redir_list: object

    def __init__(self, line_number, assignments, arguments, redir_list):
        self.line_number = line_number
        self.assignments = assignments
        self.arguments = arguments
        self.redir_list = redir_list

    def __repr__(self):
        output = "Command: {}".format(self.arguments)
        if(len(self.assignments) > 0):
            output += ", ass[{}]".format(self.assignments)
        if(len(self.redir_list) > 0):
            output += ", reds[{}]".format(self.redir_list)
        return output

    def json_serialize(self):
        json_output = make_kv(CommandNode.NodeName,
                              [self.line_number,
                               self.assignments,
                               self.arguments,
                               self.redir_list])
        return json_output

class SubshellNode(AstNode):
    NodeName = 'Subshell'
    line_number: int
    body: AstNode
    redir_list: object

    def __init__(self, line_number, body, redir_list):
        self.line_number = line_number
        self.body = body
        self.redir_list = redir_list

    def json_serialize(self):
        json_output = make_kv(SubshellNode.NodeName,
                              [self.line_number,
                               self.body,
                               self.redir_list])
        return json_output
        
class AndNode(AstNode):
    NodeName = 'And'
    left_operand: AstNode
    right_operand: AstNode

    def __init__(self, left_operand, right_operand):
        self.left_operand = left_operand
        self.right_operand = right_operand

    def __repr__(self):
        output = "{} && {}".format(self.left_operand, self.right_operand)
        return output
    
    def json_serialize(self):
        json_output = make_kv(AndNode.NodeName,
                              [self.left_operand,
                               self.right_operand])
        return json_output

class OrNode(AstNode):
    NodeName = 'Or'
    left_operand: AstNode
    right_operand: AstNode

    def __init__(self, left_operand, right_operand):
        self.left_operand = left_operand
        self.right_operand = right_operand

    def __repr__(self):
        output = "{} || {}".format(self.left_operand, self.right_operand)
        return output
    
    def json_serialize(self):
        json_output = make_kv(OrNode.NodeName,
                              [self.left_operand,
                               self.right_operand])
        return json_output
    
class SemiNode(AstNode):
    NodeName = 'Semi'
    left_operand: AstNode
    right_operand: AstNode

    def __init__(self, left_operand, right_operand):
        self.left_operand = left_operand
        self.right_operand = right_operand

    def __repr__(self):
        output = "{} ; {}".format(self.left_operand, self.right_operand)
        return output
    
    def json_serialize(self):
        json_output = make_kv(SemiNode.NodeName,
                              [self.left_operand,
                               self.right_operand])
        return json_output


class NotNode(AstNode):
    NodeName = 'Not'
    body: AstNode

    def __init__(self, body):
        self.body = body

    def json_serialize(self):
        json_output = make_kv(NotNode.NodeName,
                              self.body)
        return json_output

class RedirNode(AstNode):
    NodeName = 'Redir'
    line_number: int
    node: AstNode
    redir_list: object

    def __init__(self, line_number, node, redir_list):
        self.line_number = line_number
        self.node = node
        self.redir_list = redir_list

    def json_serialize(self):
        json_output = make_kv(RedirNode.NodeName,
                              [self.line_number,
                               self.node,
                               self.redir_list])
        return json_output

class BackgroundNode(AstNode):
    NodeName = 'Background'
    line_number: int
    node: AstNode
    redir_list: object

    def __init__(self, line_number, node, redir_list):
        self.line_number = line_number
        self.node = node
        self.redir_list = redir_list

    def json_serialize(self):
        json_output = make_kv(BackgroundNode.NodeName,
                              [self.line_number,
                               self.node,
                               self.redir_list])
        return json_output

class DefunNode(AstNode):
    NodeName = 'Defun'
    line_number: int
    name: object
    body: AstNode

    def __init__(self, line_number, name, body):
        self.line_number = line_number
        self.name = name
        self.body = body

    def json_serialize(self):
        json_output = make_kv(DefunNode.NodeName,
                              [self.line_number,
                               self.name,
                               self.body])
        return json_output

class ForNode(AstNode):
    NodeName = 'For'
    line_number: int
    argument: object
    body: AstNode
    variable: object

    def __init__(self, line_number, argument, body, variable):
        self.line_number = line_number
        self.argument = argument
        self.body = body
        self.variable = variable

    def __repr__(self):
        output = "for {} in {}; do ({})".format(self.variable, self.argument, self.body)
        return output
    
    def json_serialize(self):
        json_output = make_kv(ForNode.NodeName,
                              [self.line_number,
                               self.argument,
                               self.body,
                               self.variable])
        return json_output

class WhileNode(AstNode):
    NodeName = 'While'
    test: object
    body: AstNode

    def __init__(self, test, body):
        self.test = test
        self.body = body

    def json_serialize(self):
        json_output = make_kv(WhileNode.NodeName,
                              [self.test,
                               self.body])
        return json_output

class IfNode(AstNode):
    NodeName = 'If'
    cond: object
    then_b: AstNode
    else_b: AstNode

    def __init__(self, cond, then_b, else_b):
        self.cond = cond
        self.then_b = then_b
        self.else_b = else_b

    def json_serialize(self):
        json_output = make_kv(IfNode.NodeName,
                              [self.cond,
                               self.then_b,
                               self.else_b])
        return json_output

class CaseNode(AstNode):
    NodeName = 'Case'
    line_number: int
    argument: object
    cases: object

    def __init__(self, line_number, argument, body):
        self.line_number = line_number
        self.argument = argument
        self.body = body

    def json_serialize(self):
        json_output = make_kv(CaseNode.NodeName,
                              [self.line_number,
                               self.argument,
                               self.cases])
        return json_output
    
## This function takes an object that contains a mix of untyped and typed AstNodes (yuck) 
## and turns it into untyped json-like object. It is required atm because the infrastructure that
## we have does not translate everything to its typed form at once before compiling, and therefore
## we end up with these abomination objects.
##
## Very important TODO: 
##    We need to fix this by properly defining types (based on `compiler/parser/ast_atd.atd`)
##    and creating a bidirectional transformation from these types to the untyped json object.
##    Then we can have all ast_to_ir infrastructure work on these objects, and only in the end
##    requiring to go to the untyped form to interface with printing and parsing 
##    (which ATM does not interface with the typed form).
def ast_node_to_untyped_deep(node):
    if(isinstance(node, AstNode)):
        json_key, json_val = node.json_serialize()
        untyped_json_val = [ast_node_to_untyped_deep(obj) for obj in json_val]
        return [json_key, untyped_json_val]
    elif(isinstance(node, list)):
        return [ast_node_to_untyped_deep(obj) for obj in node]
    elif(isinstance(node, tuple)):
        return [ast_node_to_untyped_deep(obj) for obj in node]
    elif(isinstance(node, dict)):
        return {k: ast_node_to_untyped_deep(v) for k, v in node.items()}
    else:
        return node
