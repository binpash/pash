from json import JSONEncoder

from shell_ast.ast_node_c import *
from definitions.no_match_exception import *
from util import *


class AstNode:
    NodeName = 'None'

class CustomJSONEncoder(JSONEncoder):
    def default(self, obj):
        if isinstance(obj, AstNode):
            return obj.json_serialize()
        # Let the base class default method raise the TypeError
        return JSONEncoder.default(self, obj)


class PipeNode(AstNode):
    NodeName = 'Pipe'
    is_background: bool
    items: "list[AstNode]"

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
    assignments: list
    arguments: list
    redir_list: list

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
    redir_list: list

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
    redir_list: list

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
    redir_list: list

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
    test: AstNode
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
    cond: AstNode
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
    cases: list

    def __init__(self, line_number, argument, cases):
        self.line_number = line_number
        self.argument = argument
        self.cases = cases

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

def make_typed_semi_sequence(asts: "list[AstNode]") -> SemiNode:
    assert(len(asts) > 0)

    if(len(asts) == 1):
        return asts[0]
    else:
        acc = asts[-1]
        ## Remove the last ast
        iter_asts = asts[:-1]
        for ast in iter_asts[::-1]:
            acc = SemiNode(ast, acc)
        return acc
