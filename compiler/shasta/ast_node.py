import abc
from json import JSONEncoder
from shasta.print_lib import * 

class AstNode(metaclass=abc.ABCMeta):
    NodeName = 'None'

    @abc.abstractmethod
    def json(self):
        return
    
    @abc.abstractmethod
    def pretty(self) -> str:
        """
        Renders an AST back in shell syntax. 
        """
        return

class Command(AstNode):
    pass

class CustomJSONEncoder(JSONEncoder):
    def default(self, obj):
        if isinstance(obj, AstNode):
            return obj.json()
        # Let the base class default method raise the TypeError
        return JSONEncoder.default(self, obj)

class PipeNode(Command):
    NodeName = 'Pipe'
    is_background: bool
    items: "list[Command]"

    def __init__(self, is_background, items):
        self.is_background = is_background
        self.items = items

    def __repr__(self):
        if (self.is_background):
            return "Background Pipe: {}".format(self.items)    
        else:
            return "Pipe: {}".format(self.items)
        
    def json(self):
        json_output = make_kv(PipeNode.NodeName,
                              [self.is_background,
                               self.items])
        return json_output

    def pretty(self):
        bg = self.is_background
        ps = self.items
        p = intercalate(" | ", [item.pretty() for item in ps])

        if bg:
            return background(p)
        else:
            return p


class CommandNode(Command):
    NodeName = 'Command'
    line_number: int
    assignments: list
    arguments: "list[list[ArgChar]]"
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

    def json(self):
        json_output = make_kv(CommandNode.NodeName,
                              [self.line_number,
                               self.assignments,
                               self.arguments,
                               self.redir_list])
        return json_output
    
    def pretty(self):
        assigns = self.assignments
        cmds = self.arguments
        redirs = self.redir_list

        str = " ".join([assign.pretty() for assign in assigns])
        if (len(assigns) == 0) or (len(cmds) == 0):
            pass
        else:
            str += " "
        str += separated(string_of_arg, cmds) + string_of_redirs(redirs)

        return str

class SubshellNode(Command):
    NodeName = 'Subshell'
    line_number: int
    body: Command
    redir_list: list

    def __init__(self, line_number, body, redir_list):
        self.line_number = line_number
        self.body = body
        self.redir_list = redir_list

    def json(self):
        json_output = make_kv(SubshellNode.NodeName,
                              [self.line_number,
                               self.body,
                               self.redir_list])
        return json_output
    
    def pretty(self):
        return parens(self.body.pretty() + string_of_redirs(self.redir_list))
        
class AndNode(Command):
    NodeName = 'And'
    left_operand: Command
    right_operand: Command

    def __init__(self, left_operand, right_operand):
        self.left_operand = left_operand
        self.right_operand = right_operand

    def __repr__(self):
        output = "{} && {}".format(self.left_operand, self.right_operand)
        return output
    
    def json(self):
        json_output = make_kv(AndNode.NodeName,
                              [self.left_operand,
                               self.right_operand])
        return json_output
    
    def pretty(self):
        return f'{braces(self.left_operand.pretty())} && {braces(self.right_operand.pretty())}'

class OrNode(Command):
    NodeName = 'Or'
    left_operand: Command
    right_operand: Command

    def __init__(self, left_operand, right_operand):
        self.left_operand = left_operand
        self.right_operand = right_operand

    def __repr__(self):
        output = "{} || {}".format(self.left_operand, self.right_operand)
        return output
    
    def json(self):
        json_output = make_kv(OrNode.NodeName,
                              [self.left_operand,
                               self.right_operand])
        return json_output
    
    def pretty(self):
        return f'{braces(self.left_operand.pretty())} || {braces(self.right_operand.pretty())}'
    
class SemiNode(Command):
    NodeName = 'Semi'
    left_operand: Command
    right_operand: Command

    def __init__(self, left_operand, right_operand):
        self.left_operand = left_operand
        self.right_operand = right_operand

    def __repr__(self):
        output = "{} ; {}".format(self.left_operand, self.right_operand)
        return output
    
    def json(self):
        json_output = make_kv(SemiNode.NodeName,
                              [self.left_operand,
                               self.right_operand])
        return json_output
    
    def pretty(self):
        return f'{braces(self.left_operand.pretty())} \n {braces(self.right_operand.pretty())}'


class NotNode(Command):
    NodeName = 'Not'
    body: Command

    def __init__(self, body):
        self.body = body

    def json(self):
        json_output = make_kv(NotNode.NodeName,
                              self.body)
        return json_output
    
    def pretty(self):
        return f'! {braces(self.body.pretty())}'

class RedirNode(Command):
    NodeName = 'Redir'
    line_number: int
    node: Command
    redir_list: list

    def __init__(self, line_number, node, redir_list):
        self.line_number = line_number
        self.node = node
        self.redir_list = redir_list

    def json(self):
        json_output = make_kv(RedirNode.NodeName,
                              [self.line_number,
                               self.node,
                               self.redir_list])
        return json_output
    
    def pretty(self):
        return self.node.pretty() + string_of_redirs(self.redir_list)

class BackgroundNode(Command):
    NodeName = 'Background'
    line_number: int
    node: Command
    redir_list: list

    def __init__(self, line_number, node, redir_list):
        self.line_number = line_number
        self.node = node
        self.redir_list = redir_list

    def json(self):
        json_output = make_kv(BackgroundNode.NodeName,
                              [self.line_number,
                               self.node,
                               self.redir_list])
        return json_output

    def pretty(self):
        return background(self.node.pretty() + string_of_redirs(self.redir_list))

class DefunNode(Command):
    NodeName = 'Defun'
    line_number: int
    name: object
    body: Command

    def __init__(self, line_number, name, body):
        self.line_number = line_number
        self.name = name
        self.body = body

    def json(self):
        json_output = make_kv(DefunNode.NodeName,
                              [self.line_number,
                               self.name,
                               self.body])
        return json_output
    
    def pretty(self):
        name = self.name
        body = self.body
        return name + "() {\n" + body.pretty() + "\n}"


class ForNode(Command):
    NodeName = 'For'
    line_number: int
    argument: "list[list[ArgChar]]"
    body: Command
    variable: object

    def __init__(self, line_number, argument, body, variable):
        self.line_number = line_number
        self.argument = argument
        self.body = body
        self.variable = variable

    def __repr__(self):
        output = "for {} in {}; do ({})".format(self.variable, self.argument, self.body)
        return output
    
    def json(self):
        json_output = make_kv(ForNode.NodeName,
                              [self.line_number,
                               self.argument,
                               self.body,
                               self.variable])
        return json_output

    def pretty(self):
        a = self.argument
        body = self.body
        var = self.variable
        return f'for {var} in {separated(string_of_arg, a)}; do {body.pretty()}; done'

class WhileNode(Command):
    NodeName = 'While'
    test: Command
    body: Command

    def __init__(self, test, body):
        self.test = test
        self.body = body

    def json(self):
        json_output = make_kv(WhileNode.NodeName,
                              [self.test,
                               self.body])
        return json_output

    def pretty(self):
        first = self.test
        b = self.body
        
        if isinstance(first, NotNode):
            t = first.body
            return f'until {t.pretty()}; do {b.pretty()}; done '
        else:
            t = first
            return f'while {t.pretty()}; do {b.pretty()}; done '

class IfNode(Command):
    NodeName = 'If'
    cond: Command
    then_b: Command
    else_b: Command

    def __init__(self, cond, then_b, else_b):
        self.cond = cond
        self.then_b = then_b
        self.else_b = else_b

    def json(self):
        json_output = make_kv(IfNode.NodeName,
                              [self.cond,
                               self.then_b,
                               self.else_b])
        return json_output

    def pretty(self):
        c = self.cond
        t = self.then_b
        e = self.else_b
        str1 = f'if {c.pretty()}; then {t.pretty()}'

        if is_empty_cmd(e):
            str1 += "; fi"
        elif isinstance(e, IfNode):
            str1 += "; el" + e.pretty()
        else:
            str1 += f'; else {e.pretty()}; fi'

        return str1

class CaseNode(Command):
    NodeName = 'Case'
    line_number: int
    argument: "list[ArgChar]"
    cases: list

    def __init__(self, line_number, argument, cases):
        self.line_number = line_number
        self.argument = argument
        self.cases = cases

    def json(self):
        json_output = make_kv(CaseNode.NodeName,
                              [self.line_number,
                               self.argument,
                               self.cases])
        return json_output
    
    def pretty(self):
        a = self.argument
        cs = self.cases
        return f'case {string_of_arg(a)} in {separated(string_of_case, cs)} esac'


class ArgChar(AstNode):
    ## This method formats an arg_char to a string to
    ## the best of its ability
    def format(self) -> str:
        raise NotImplementedError

class CArgChar(ArgChar):
    NodeName = 'C'
    char: int

    def __init__(self, char: int):
        self.char = char

    def __repr__(self):
        return self.format()
    
    def format(self) -> str:
        return str(chr(self.char))

    def json(self):
        json_output = make_kv(CArgChar.NodeName,
                              self.char)
        return json_output
    
    def pretty(self, quote_mode=UNQUOTED):
        if quote_mode==QUOTED and chr(self.char) == '"':
            return '\\"'
        else:
            return chr(self.char)

class EArgChar(ArgChar):
    NodeName = 'E'
    char: int

    def __init__(self, char: int):
        self.char = char

    ## TODO: Implement
    def __repr__(self):
        return f'\\{chr(self.char)}'

    def format(self) -> str:
        ## TODO: This is not right. I think the main reason for the
        ## problems is the differences between bash and the posix
        ## standard.
        non_escape_chars = [92, # \
                            61, # =
                            91, # [
                            93, # ]
                            45, # -
                            58, # :
                            126,# ~
                            42] # *
        if(self.char in non_escape_chars):
            return '{}'.format(chr(self.char))
        else:
            return '\{}'.format(chr(self.char))

    def json(self):
        json_output = make_kv(EArgChar.NodeName,
                              self.char)
        return json_output
    
    def pretty(self, quote_mode=UNQUOTED):
        param = self.char
        char = chr(param)

        ## MMG 2021-09-20 It might be safe to move everything except for " in the second list, but no need to do it if the tests pass
        ## '!' dropped for bash non-interactive bash compatibility
        ## Chars to escape unconditionally
        chars_to_escape = ["'", '"', '`', '(', ')', '{', '}', '$', '&', '|', ';']
        ## Chars to escape only when not quoted
        chars_to_escape_when_no_quotes = ['*', '?', '[', ']', '#', '<', '>', '~', ' ']
        if char in chars_to_escape:
            return '\\' + char
        elif char in chars_to_escape_when_no_quotes and quote_mode==UNQUOTED:
            return '\\' + char
        else:
            return escaped(param)


class TArgChar(ArgChar):
    NodeName = 'T'
    string: str

    def __init__(self, string: str):
        self.string = string

    ## TODO: Implement
    # def __repr__(self):
    #     return f''

    def json(self):
        json_output = make_kv(TArgChar.NodeName,
                              self.string)
        return json_output

    def pretty(self, quote_mode=UNQUOTED):
        param = self.string
        ## TODO: Debug this
        if param == "None":
            return "~"
        elif len(param) == 2:
            if param[0] == "Some":
                (_, u) = param

                return "~" + u
            else:
                assert(False)
        else:
            print ("Unexpected param for T: %s" % param)

class AArgChar(ArgChar):
    NodeName = 'A'
    arg: "list[ArgChar]"

    def __init__(self, arg: "list[ArgChar]"):
        self.arg = arg

    ## TODO: Implement
    # def __repr__(self):
    #     return f''

    def json(self):
        json_output = make_kv(AArgChar.NodeName,
                              self.arg)
        return json_output
    
    def pretty(self, quote_mode=UNQUOTED):
        param = self.arg
        return f'$(({string_of_arg(param, quote_mode)}))'

class VArgChar(ArgChar):
    NodeName = 'V'
    fmt: object
    null: bool
    var: str
    arg: "list[ArgChar]"

    def __init__(self, fmt, null: bool, var: str, arg: "list[ArgChar]"):
        self.fmt = fmt
        self.null = null
        self.var = var
        self.arg = arg

    def __repr__(self):
        return f'V({self.fmt},{self.null},{self.var},{self.arg})'

    def format(self) -> str:
        return '${{{}}}'.format(self.var)

    def json(self):
        json_output = make_kv(VArgChar.NodeName,
                              [self.fmt,
                               self.null,
                               self.var,
                               self.arg])
        return json_output
    
    def pretty(self, quote_mode=UNQUOTED):
        vt = self.fmt
        nul = self.null
        name = self.var
        a = self.arg
        if vt == "Length":
            return "${#" + name + "}"
        else:
            stri = "${" + name

            # Depending on who generated the JSON, nul may be
            # a string or a boolean! In Python, non-empty strings
            # to True.
            if (str(nul).lower() == "true"):
                stri += ":"
            elif (str (nul).lower() == "false"):
                pass
            else:
                assert(False)

            stri += string_of_var_type(vt) + string_of_arg(a, quote_mode) + "}"

            return stri


class QArgChar(ArgChar):
    NodeName = 'Q'
    arg: "list[ArgChar]"

    def __init__(self, arg: "list[ArgChar]"):
        self.arg = arg

    def __repr__(self):
        return f'Q({self.arg})'
    
    def format(self) -> str:
        chars = [arg_char.format() for arg_char in self.arg]
        joined_chars = "".join(chars)
        return '"{}"'.format(joined_chars)

    def json(self):
        json_output = make_kv(QArgChar.NodeName,
                              self.arg)
        return json_output

    def pretty(self, quote_mode=UNQUOTED):
        param = self.arg
        return "\"" + string_of_arg(param, quote_mode=QUOTED) + "\""


class BArgChar(ArgChar):
    NodeName = 'B'
    node: Command

    def __init__(self, node: Command):
        self.node = node

    ## TODO: Implement
    # def __repr__(self):
    #     return f''

    def format(self) -> str:
        return '$({})'.format(self.node)

    def json(self):
        json_output = make_kv(BArgChar.NodeName,
                              self.node)
        return json_output
    
    def pretty(self, quote_mode=UNQUOTED):
        param = self.node
        return "$(" + param.pretty() + ")"

class AssignNode(AstNode):
    var: str
    val: "list[ArgChar]"

    def __init__(self, var: str, val):
        self.var = var
        self.val = val

    # TODO: Implement
    def __repr__(self):
        return f'{self.var}={self.val}'

    def json(self):
        json_output = [self.var, self.val]
        return json_output
    
    def pretty(self):
        return f'{self.var}={string_of_arg(self.val)}'

    
class RedirectionNode(AstNode):
    redir_type: str
    fd: int
    arg: "list[ArgChar]"
    pass

class FileRedirNode(RedirectionNode):
    NodeName = "File"
    redir_type: str
    fd: int
    arg: "list[ArgChar]"

    def __init__(self, redir_type, fd, arg):
        self.redir_type = redir_type
        self.fd = fd
        self.arg = arg

    # TODO: Implement
    # def __repr__(self):
    #     return f''

    def json(self):
        json_output = make_kv(FileRedirNode.NodeName,
                              [self.redir_type,
                               self.fd,
                               self.arg])
        return json_output
    
    def pretty(self):
        subtype = self.redir_type
        fd = self.fd
        a = self.arg
        if subtype == "To":
            return show_unless(1, fd) + ">" + string_of_arg(a)
        elif subtype == "Clobber":
            return show_unless(1, fd) + ">|" + string_of_arg(a)
        elif subtype == "From":
            return show_unless(0, fd) + "<" + string_of_arg(a)
        elif subtype == "FromTo":
            return show_unless(0, fd) + "<>" + string_of_arg(a)
        elif subtype == "Append":
            return show_unless(1, fd) + ">>" + string_of_arg(a)
        assert(False)


class DupRedirNode(RedirectionNode):
    NodeName = "Dup"
    dup_type: str
    fd: int
    arg: "list[ArgChar]"

    def __init__(self, dup_type, fd, arg):
        self.dup_type = dup_type
        self.fd = fd
        self.arg = arg

    # TODO: Implement
    # def __repr__(self):
    #     return f''

    def json(self):
        json_output = make_kv(DupRedirNode.NodeName,
                              [self.dup_type,
                               self.fd,
                               self.arg])
        return json_output
    
    def pretty(self):
        subtype = self.dup_type
        fd = self.fd
        tgt = self.arg
        if subtype == "ToFD":
            return show_unless(1, fd) + ">&" + string_of_arg(tgt)
        elif subtype == "FromFD":
            return show_unless(0, fd) + "<&" + string_of_arg(tgt)
        assert(False)
    
class HeredocRedirNode(RedirectionNode):
    NodeName = "Heredoc"
    heredoc_type: str
    fd: int
    arg: "list[ArgChar]"

    def __init__(self, heredoc_type, fd, arg):
        self.heredoc_type = heredoc_type
        self.fd = fd
        self.arg = arg

    # TODO: Implement
    # def __repr__(self):
    #     return f''

    def json(self):
        json_output = make_kv(HeredocRedirNode.NodeName,
                              [self.heredoc_type,
                               self.fd,
                               self.arg])
        return json_output
    
    def pretty(self):
        t = self.heredoc_type
        fd = self.fd
        a = self.arg
        heredoc = string_of_arg(a, quote_mode=HEREDOC)
        marker = fresh_marker0(heredoc)

        stri = show_unless(0, fd) + "<<"
        if t == "XHere":
            stri += marker
        else:
            stri += "'" + marker + "'"

        stri += "\n" + heredoc + marker + "\n"

        return stri

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
        json_key, json_val = node.json()
        return [json_key, ast_node_to_untyped_deep(json_val)]
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

## This has to be here and not in print_lib to avoid circular dependencies
def string_of_arg(args, quote_mode=UNQUOTED):
    i = 0
    text = []
    while i < len(args):
        c = args[i].pretty(quote_mode=quote_mode)
        if c == "$" and (i+1 < len(args)) and isinstance(args[i+1],EArgChar):
            c = "\\$"
        text.append(c)

        i = i+1
    
    text = "".join(text)

    return text

def string_of_case(c):
    pats = map(string_of_arg, c["cpattern"])

    return f'{intercalate("|", pats)}) {c["cbody"].pretty()};;'


def is_empty_cmd(e: Command):
    return isinstance(e, CommandNode) \
        and e.line_number == -1 \
        and len(e.assignments) == 0 \
        and len(e.arguments) == 0 \
        and len(e.redir_list) == 0


## Implements a pattern-matching style traversal over the AST
def ast_match(ast_node, cases, *args):
    return cases[type(ast_node).NodeName](*args)(ast_node)

## Util function
def make_kv(key, val):
    return [key, val]
