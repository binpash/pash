
from shell_ast.ast_node import *

## TODO: Recurse in necessary subfields and iterate deep
def to_ast_node(obj) -> AstNode:
    k, v = obj
    if k == PipeNode.NodeName:
        node = PipeNode(is_background=v[0],
                        items=v[1])
    elif k == CommandNode.NodeName:
        node = CommandNode(line_number=v[0],
                           assignments=v[1],
                           arguments=v[2],
                           redir_list=v[3])
    elif k == SubshellNode.NodeName:
        node = SubshellNode(line_number=v[0],
                            body=v[1],
                            redir_list=v[2])
    elif k == AndNode.NodeName:
        node = AndNode(left_operand=v[0],
                        right_operand=v[1])
    elif k == OrNode.NodeName:
        node = OrNode(left_operand=v[0],
                        right_operand=v[1])
    elif k == SemiNode.NodeName:
        node = SemiNode(left_operand=v[0],
                        right_operand=v[1])
    elif k == NotNode.NodeName:
        node = NotNode(body=v)
    elif k == RedirNode.NodeName:
        node = RedirNode(line_number=v[0],
                         node=v[1],
                         redir_list=v[2])
    elif k == BackgroundNode.NodeName:
        node = BackgroundNode(line_number=v[0],
                              node=v[1],
                              redir_list=v[2])
    elif k == DefunNode.NodeName:
        node = DefunNode(line_number=v[0],
                         name=v[1],
                         body=v[2])
    elif k == ForNode.NodeName:
        node = ForNode(line_number=v[0],
                       argument=v[1],
                       body=v[2],
                       variable=v[3])
    elif k == WhileNode.NodeName:
        node = WhileNode(test=v[0],
                         body=v[1])
    elif k == IfNode.NodeName:
        node = IfNode(cond=v[0],
                      then_b=v[1],
                      else_b=v[2])
    elif k == CaseNode.NodeName:
        node = CaseNode(line_number=v[0],
                        argument=v[1],
                        cases=v[2])
    else:
        raise ValueError()
    return node
