"""
Generic AST walker for shell AST nodes.

Provides a reusable tree traversal with visitor and replacement callbacks.
"""

from shasta.ast_node import (
    # Base classes
    AstNode,
    ArgChar,
    # ArgChar types
    CArgChar,
    EArgChar,
    TArgChar,
    AArgChar,
    VArgChar,
    QArgChar,
    BArgChar,
    # Command nodes
    PipeNode,
    CommandNode,
    SubshellNode,
    AndNode,
    OrNode,
    SemiNode,
    NotNode,
    RedirNode,
    BackgroundNode,
    DefunNode,
    ForNode,
    WhileNode,
    IfNode,
    CaseNode,
    # Redirection nodes
    FileRedirNode,
    DupRedirNode,
    HeredocRedirNode,
    # Assignment
    AssignNode,
    # Bash-specific nodes
    SelectNode,
    ArithNode,
    CondNode,
    ArithForNode,
    CoprocNode,
    TimeNode,
    SingleArgRedirNode,
    GroupNode,
)


def walk_ast_node(node, visit=None, replace=None):
    """
    Generic preorder visitor/transformer for shell AST nodes.

    This function traverses the AST in preorder (parent before children),
    optionally calling visitor and replacement callbacks at each node.

    :param node: The AST node to walk
    :param visit: Optional visitor callback called for side effects (return value ignored)
    :param replace: Optional replacement callback; if it returns non-None, that value
                   replaces the node and recursion stops for that subtree
    :returns: The (potentially transformed) node
    """
    if visit:
        visit(node)
    if replace:
        replaced = replace(node)
        if replaced is not None:
            return replaced

    def walk(n):
        return walk_ast_node(n, visit=visit, replace=replace)

    def walk_fd(fd):
        """Walk file descriptor which can be ('var', argchars) or ('fixed', int)."""
        match fd:
            case ("var", argchars):
                return ("var", walk(argchars))
            case _:
                return fd

    match node:
        # Handle collections
        case list():
            return [walk(n) for n in node]
        case tuple():
            return tuple(walk(n) for n in node)

        # ArgChar types
        case CArgChar() | EArgChar() | TArgChar():
            # Leaf nodes with no children to walk
            return node
        case AArgChar():
            return AArgChar(
                arg=walk(node.arg),
            )
        case VArgChar():
            return VArgChar(
                fmt=node.fmt,
                null=node.null,
                var=node.var,
                arg=walk(node.arg),
            )
        case QArgChar():
            return QArgChar(
                arg=walk(node.arg),
            )
        case BArgChar():
            return BArgChar(
                node=walk(node.node),
            )

        # Command nodes
        case PipeNode():
            return PipeNode(
                is_background=node.is_background,
                items=[walk(item) for item in node.items],
            )
        case CommandNode():
            return CommandNode(
                line_number=node.line_number,
                assignments=[walk(a) for a in node.assignments],
                arguments=[walk(arg) for arg in node.arguments],
                redir_list=[walk(r) for r in node.redir_list],
            )
        case SubshellNode():
            return SubshellNode(
                line_number=node.line_number,
                body=walk(node.body),
                redir_list=[walk(r) for r in node.redir_list],
            )
        case AndNode():
            return AndNode(
                left_operand=walk(node.left_operand),
                right_operand=walk(node.right_operand),
                no_braces=node.no_braces,
            )
        case OrNode():
            return OrNode(
                left_operand=walk(node.left_operand),
                right_operand=walk(node.right_operand),
                no_braces=node.no_braces,
            )
        case SemiNode():
            return SemiNode(
                left_operand=walk(node.left_operand),
                right_operand=walk(node.right_operand),
                semicolon=node.semicolon,
            )
        case NotNode():
            return NotNode(
                body=walk(node.body),
                no_braces=node.no_braces,
            )
        case RedirNode():
            return RedirNode(
                line_number=node.line_number,
                node=walk(node.node),
                redir_list=[walk(r) for r in node.redir_list],
            )
        case BackgroundNode():
            return BackgroundNode(
                line_number=node.line_number,
                node=walk(node.node),
                redir_list=[walk(r) for r in node.redir_list],
                after_ampersand=walk(node.after_ampersand) if node.after_ampersand else None,
                no_braces=node.no_braces,
            )
        case DefunNode():
            return DefunNode(
                line_number=node.line_number,
                name=walk(node.name),
                body=walk(node.body),
                bash_mode=node.bash_mode,
            )
        case ForNode():
            return ForNode(
                line_number=node.line_number,
                argument=[walk(arg) for arg in node.argument],
                body=walk(node.body),
                variable=walk(node.variable),
            )
        case WhileNode():
            return WhileNode(
                test=walk(node.test),
                body=walk(node.body),
            )
        case IfNode():
            return IfNode(
                cond=walk(node.cond),
                then_b=walk(node.then_b),
                else_b=walk(node.else_b) if node.else_b else None,
            )
        case CaseNode():
            updated_cases = []
            for case in node.cases:
                new_case = dict(case)
                if "cpattern" in case:
                    new_case["cpattern"] = [walk(p) for p in case["cpattern"]]
                if case.get("cbody"):
                    new_case["cbody"] = walk(case["cbody"])
                updated_cases.append(new_case)
            return CaseNode(
                line_number=node.line_number,
                argument=walk(node.argument),
                cases=updated_cases,
            )

        # Assignment node
        case AssignNode():
            return AssignNode(
                var=node.var,
                val=walk(node.val),
            )

        # Redirection nodes
        case FileRedirNode():
            return FileRedirNode(
                redir_type=node.redir_type,
                fd=walk_fd(node.fd),
                arg=walk(node.arg) if node.arg else None,
            )
        case DupRedirNode():
            return DupRedirNode(
                dup_type=node.dup_type,
                fd=walk_fd(node.fd),
                arg=walk_fd(node.arg),
                move=node.move,
            )
        case HeredocRedirNode():
            return HeredocRedirNode(
                heredoc_type=node.heredoc_type,
                fd=walk_fd(node.fd),
                arg=walk(node.arg),
                kill_leading=node.kill_leading,
                eof=node.eof,
            )
        case SingleArgRedirNode():
            return SingleArgRedirNode(
                redir_type=node.redir_type,
                fd=walk_fd(node.fd),
            )

        # Bash-specific command nodes
        case SelectNode():
            return SelectNode(
                line_number=node.line_number,
                variable=walk(node.variable),
                body=walk(node.body),
                map_list=[walk(m) for m in node.map_list],
            )
        case ArithNode():
            return ArithNode(
                line_number=node.line_number,
                body=[walk(b) for b in node.body],
            )
        case CondNode():
            return CondNode(
                line_number=node.line_number,
                cond_type=node.cond_type,
                op=walk(node.op) if node.op else None,
                left=walk(node.left) if node.left else None,
                right=walk(node.right) if node.right else None,
                invert_return=node.invert_return,
            )
        case ArithForNode():
            return ArithForNode(
                line_number=node.line_number,
                init=walk(node.init),
                cond=walk(node.cond),
                step=walk(node.step),
                action=walk(node.action),
            )
        case CoprocNode():
            return CoprocNode(
                name=walk(node.name),
                body=walk(node.body),
            )
        case TimeNode():
            return TimeNode(
                time_posix=node.time_posix,
                command=walk(node.command),
            )
        case GroupNode():
            return GroupNode(
                body=walk(node.body),
            )

        # Default: return node unchanged (handles primitives, None, etc.)
        case _:
            return node
