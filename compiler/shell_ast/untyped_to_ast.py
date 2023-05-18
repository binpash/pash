
from shell_ast.ast_node import *

## TODO: Recurse in necessary subfields and iterate deep
def to_ast_node(obj) -> AstNode:
    k, v = obj
    if k == PipeNode.NodeName:
        node = PipeNode(is_background=v[0],
                        items=to_ast_nodes(v[1]))
    elif k == CommandNode.NodeName:
        node = CommandNode(line_number=v[0],
                           assignments=to_assigns(v[1]),
                           arguments=to_args(v[2]),
                           redir_list=to_redirs(v[3]))
    elif k == SubshellNode.NodeName:
        node = SubshellNode(line_number=v[0],
                            body=to_ast_node(v[1]),
                            redir_list=to_redirs(v[2]))
    elif k == AndNode.NodeName:
        node = AndNode(left_operand=to_ast_node(v[0]),
                       right_operand=to_ast_node(v[1]))
    elif k == OrNode.NodeName:
        node = OrNode(left_operand=to_ast_node(v[0]),
                      right_operand=to_ast_node(v[1]))
    elif k == SemiNode.NodeName:
        node = SemiNode(left_operand=to_ast_node(v[0]),
                        right_operand=to_ast_node(v[1]))
    elif k == NotNode.NodeName:
        node = NotNode(body=to_ast_node(v))
    elif k == RedirNode.NodeName:
        node = RedirNode(line_number=v[0],
                         node=to_ast_node(v[1]),
                         redir_list=to_redirs(v[2]))
    elif k == BackgroundNode.NodeName:
        node = BackgroundNode(line_number=v[0],
                              node=to_ast_node(v[1]),
                              redir_list=to_redirs(v[2]))
    elif k == DefunNode.NodeName:
        node = DefunNode(line_number=v[0],
                         name=v[1],
                         body=to_ast_node(v[2]))
    elif k == ForNode.NodeName:
        node = ForNode(line_number=v[0],
                       argument=to_args(v[1]),
                       body=to_ast_node(v[2]),
                       variable=v[3])
    elif k == WhileNode.NodeName:
        node = WhileNode(test=to_ast_node(v[0]),
                         body=to_ast_node(v[1]))
    elif k == IfNode.NodeName:
        node = IfNode(cond=to_ast_node(v[0]),
                      then_b=to_ast_node(v[1]),
                      else_b=to_ast_node(v[2]))
    elif k == CaseNode.NodeName:
        node = CaseNode(line_number=v[0],
                        argument=to_arg(v[1]),
                        cases=to_case_list(v[2]))
    else:
        raise ValueError()
    return node

def to_ast_nodes(node_list):
    new_node_list = [to_ast_node(ast_node) for ast_node in node_list]
    return new_node_list

def to_assigns(assignments):
    new_assignments = []
    for name, val in assignments:
        new_assignments.append(make_kv(name, to_arg(val)))
    return new_assignments

def to_redirs(redir_list):
    new_redir_list = [to_redir(redir) for redir in redir_list]
    return new_redir_list

def to_redir(redir):
    k, v = redir
    if k == "File":
        return make_kv(k, 
                       [v[0],
                        v[1],
                        to_arg(v[2])])
    elif k == "Dup":
        return make_kv(k, 
                       [v[0],
                        v[1],
                        to_arg(v[2])])
    elif k == "Heredoc":
        return make_kv(k, 
                       [v[0],
                        v[1],
                        to_arg(v[2])])
    assert(False)

def to_args(arg_list):
    new_arg_list = [to_arg(arg) for arg in arg_list]
    return new_arg_list

## TODO: Make this return an Arg (since we have that class already)
def to_arg(arg_char_list):
    new_arg_char_list = [to_arg_char(arg_char) for arg_char in arg_char_list]
    return new_arg_char_list

def to_arg_char(arg_char):
    k, v = arg_char
    if k == "C":
        return arg_char
    elif k == "E":
        return arg_char
    elif k == "T":
        return arg_char
    elif k == "A":
        return make_kv(k, to_arg(v))
    elif k == "V":
        return make_kv(k,
                       [v[0],
                        v[1],
                        v[2],
                        to_arg(v[3])])
    elif k == "Q":
        return make_kv(k, to_arg(v))
    elif k == "B":
        return make_kv(k, to_ast_node(v))
    assert(False)

def to_case_list(case_list):
    new_case_list = [(to_args(case[0]), to_ast_node(case[1])) 
                     for case in case_list]
    return new_case_list
