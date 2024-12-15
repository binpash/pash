from __future__ import annotations

from shasta.ast_node import *
import config
from util import *
import pexpect

VAR_FILE_PATH = None


def expand_using_bash(ast: list[AstNode]) -> list[AstNode]:
    global VAR_FILE_PATH
    VAR_FILE_PATH = ptempfile()
    log("VAR_FILE_PATH:", VAR_FILE_PATH)
    with open(VAR_FILE_PATH, "w") as f:
        f.write("")

    nodes = []
    for node in ast:
        nodes.append(compile_node(node))

    os.remove(VAR_FILE_PATH)

    return nodes


def init_bash_mirror_subprocess() -> pexpect.spawn:
    ## Spawn a bash process to ask it for expansions
    p = pexpect.spawn("/usr/bin/env", ["bash", "-i"], encoding="utf-8", echo=False)
    ## If we are in debug mode also log the bash's output
    if config.pash_args.debug >= 1:
        _, file_to_save_output = ptempfile()
        log("bash mirror log saved in:", file_to_save_output)
        fout = open(file_to_save_output, "w", encoding="utf-8")
        p.logfile = fout
    return p


def query_expand_bash_mirror(string: bytes, bash_mirror):
    command = f"echo -n {string}"
    return sync_run_line_command_mirror(command, bash_mirror)


def sync_run_line_command_mirror(command: bytes, bash_mirror: pexpect.spawn):
    bash_command = command
    log("Executing bash command in mirror:", bash_command)

    # Note: this will eventually need to be changed to support non-utf8 characters
    bash_mirror.sendline(str(bash_command))

    data = wait_bash_mirror(bash_mirror)
    log("mirror done!")

    return data


def wait_bash_mirror(bash_mirror: pexpect.spawn):
    r = bash_mirror.expect(r"EXPECT\$ ")
    assert r == 0
    output: bytes = bash_mirror.before

    ## I am not sure why, but \r s are added before \n s
    output = output.replace("\r\n", "\n")

    log("Before the prompt!")
    log(output)
    return output


def update_bash_mirror_vars(bash_mirror):
    assert VAR_FILE_PATH != "" and not VAR_FILE_PATH is None

    bash_mirror.sendline(f'PS1="EXPECT\$ "')
    wait_bash_mirror(bash_mirror)
    log("PS1 set!")

    ## TODO: There is unnecessary write/read to this var file now.
    bash_mirror.sendline(f"source {VAR_FILE_PATH}")
    log("sent source to mirror")
    wait_bash_mirror(bash_mirror)
    log("mirror done!")


# ---------------------


def compile_node(ast_object: AstNode) -> AstNode:
    node_name = ast_object.NodeName

    if node_name == "Pipe":
        return compile_node_pipe(ast_object)
    elif node_name == "Command":
        return compile_node_command(ast_object)
    elif node_name == "Subshell":
        return compile_node_subshell(ast_object)
    elif node_name == "And":
        return compile_node_and(ast_object)
    elif node_name == "Or":
        return compile_node_or(ast_object)
    elif node_name == "Semi":
        return compile_node_semi(ast_object)
    elif node_name == "Not":
        return compile_node_not(ast_object)
    elif node_name == "Redir":
        return compile_node_redir(ast_object)
    elif node_name == "Background":
        return compile_node_background(ast_object)
    elif node_name == "Defun":
        return compile_node_defun(ast_object)
    elif node_name == "For":
        return compile_node_for(ast_object)
    elif node_name == "While":
        return compile_node_while(ast_object)
    elif node_name == "If":
        return compile_node_if(ast_object)
    elif node_name == "Case":
        return compile_node_case(ast_object)
    elif node_name == "Select":
        return compile_node_select(ast_object)
    elif node_name == "Arith":
        return compile_node_arith(ast_object)
    elif node_name == "Cond":
        return compile_node_cond(ast_object)
    elif node_name == "ArithFor":
        return compile_node_arith_for(ast_object)
    elif node_name == "Coproc":
        return compile_node_coproc(ast_object)
    elif node_name == "Time":
        return compile_node_time(ast_object)
    elif node_name == "Group":
        return compile_node_group(ast_object)
    else:
        log(f"Unknown node: {node_name}")
        print(f"Unknown node: {node_name}")
        raise NotImplementedError()


def compile_node_pipe(ast_node: PipeNode):
    ast_node.items = [compile_node(item) for item in ast_node.items]
    return ast_node


def compile_node_command(ast_node: CommandNode):
    ast_node.assignments = compile_command_assignments(ast_node.assignments)
    ast_node.arguments = compile_command_arguments(ast_node.arguments)
    ast_node.redir_list = compile_redirections(ast_node.redir_list)
    return ast_node


def compile_node_subshell(ast_node: SubshellNode):
    ast_node.body = compile_node(ast_node.body)
    return ast_node


def compile_node_and(ast_node: AndNode):
    ast_node.left_operand = compile_node(ast_node.left_operand)
    ast_node.right_operand = compile_node(ast_node.right_operand)
    return ast_node


def compile_node_or(ast_node: OrNode):
    ast_node.left_operand = compile_node(ast_node.left_operand)
    ast_node.right_operand = compile_node(ast_node.right_operand)
    return ast_node


def compile_node_semi(ast_node: SemiNode):
    ast_node.left_operand = compile_node(ast_node.left_operand)
    ast_node.right_operand = compile_node(ast_node.right_operand)
    return ast_node


def compile_node_not(ast_node: NotNode):
    ast_node.body = compile_node(ast_node.body)
    return ast_node


def compile_node_redir(ast_node: RedirNode):
    ast_node.node = compile_node(ast_node.node)
    ast_node.redir_list = compile_redirections(ast_node.redir_list)
    return ast_node


def compile_node_background(ast_node: BackgroundNode):
    ast_node.node = compile_node(ast_node.node)
    ast_node.redir_list = compile_redirections(ast_node.redir_list)
    return ast_node


def compile_node_defun(ast_node: DefunNode):
    ast_node.name = compile_command_argument(ast_node.name)
    ast_node.body = compile_node(ast_node.body)
    return ast_node


def compile_node_for(ast_node: ForNode):
    ast_node.variable = compile_command_argument(ast_node.variable)
    ast_node.argument = compile_command_arguments(ast_node.argument)
    # we can't do this because of argument expansion
    # ast_node.body = compile_node(ast_node.body)
    return ast_node


def compile_node_while(ast_node: WhileNode):
    ast_node.test = compile_command_argument(ast_node.test)
    ast_node.body = compile_node(ast_node.body)
    return ast_node


def compile_node_if(ast_node: IfNode):
    ast_node.cond = compile_command_argument(ast_node.cond)
    ast_node.then_b = compile_node(ast_node.then_b)
    ast_node.else_b = compile_node(ast_node.else_b) if ast_node.else_b else None
    return ast_node


def compile_node_case(ast_node: CaseNode):
    ast_node.argument = compile_command_argument(ast_node.argument)
    ast_node.cases = compile_command_cases(ast_node.cases)
    return ast_node


def compile_node_select(ast_node: SelectNode):
    ast_node.variable = compile_command_argument(ast_node.variable)
    ast_node.body = compile_node(ast_node.body)
    ast_node.map_list = compile_command_arguments(ast_node.map_list)
    return ast_node


def compile_node_arith(ast_node: ArithNode):
    ast_node.body = compile_command_arguments(ast_node.body)
    return ast_node


def compile_node_cond(ast_node: CondNode):
    ast_node.op = compile_command_argument(ast_node.op) if ast_node.op else None
    ast_node.left = compile_command_argument(ast_node.left) if ast_node.left else None
    ast_node.right = (
        compile_command_argument(ast_node.right) if ast_node.right else None
    )
    return ast_node


def compile_node_arith_for(ast_node: ArithForNode):
    ast_node.init = compile_command_arguments(ast_node.init)
    ast_node.cond = compile_command_arguments(ast_node.cond)
    ast_node.step = compile_command_arguments(ast_node.step)
    ast_node.action = compile_node(ast_node.action)
    return ast_node


def compile_node_coproc(ast_node: CoprocNode):
    ast_node.name = compile_command_argument(ast_node.name)
    ast_node.body = compile_node(ast_node.body)
    return ast_node


def compile_node_time(ast_node: TimeNode):
    ast_node.command = compile_node(ast_node.command)
    return ast_node


def compile_node_group(ast_node: GroupNode):
    ast_node.body = compile_node(ast_node.body)
    ast_node.redirections = compile_redirections(ast_node.redirections)
    return ast_node


def compile_command_arguments(arguments: List[List[ArgChar]]) -> List[List[ArgChar]]:
    return [compile_command_argument(argument) for argument in arguments]


def compile_command_argument(argument: List[CArgChar]) -> List[CArgChar]:
    bash_mirror = init_bash_mirror_subprocess()
    update_bash_mirror_vars(bash_mirror)
    expanded = query_expand_bash_mirror(
        "".join([chr(c.char) for c in argument]), bash_mirror
    )
    bash_mirror.close()
    node = bytes_to_arg_char_list(expanded)
    return node


def bytes_to_arg_char_list(data: str) -> List[ArgChar]:
    return [CArgChar(ord(c)) for c in data]


def compile_command_assignments(assignments: List[AssignNode]) -> List[AssignNode]:
    return [compile_command_assignment(assignment) for assignment in assignments]


def compile_command_assignment(assignment: AssignNode) -> AssignNode:
    assignment.val = compile_command_argument(assignment.val)
    with open(VAR_FILE_PATH, "a") as f:
        f.write(f"{assignment.var}={assignment.val}\n")
    return assignment


def compile_redirections(redir_list: List[RedirectionNode]) -> List[RedirectionNode]:
    return [compile_redirection(redir) for redir in redir_list]


def compile_redirection(redir: RedirectionNode) -> RedirectionNode:
    type = redir.redir_type
    if type == "File":
        return compile_redirection_file(redir)
    elif type == "Dup":
        return compile_redirection_dup(redir)
    elif type == "Here":
        return compile_redirection_here(redir)
    elif type == "SingleArg":
        return compile_redirection_single_arg(redir)
    else:
        log(f"Unknown redirection type: {type}")
        raise NotImplementedError()


def compile_command_cases(cases: list[dict]) -> list[dict]:
    return [compile_command_case(case) for case in cases]


def compile_command_case(case: dict) -> dict:
    case["pattern"] = compile_command_argument(case["pattern"])
    case["body"] = compile_node(case["body"])
    return case


def compile_redirection_file(redir: FileRedirNode) -> FileRedirNode:
    if redir.fd[0] == "var":
        redir.fd[1] = compile_command_argument(redir.fd[1])
    redir.arg = compile_command_argument(redir.arg)
    return redir


def compile_redirection_dup(redir: DupRedirNode) -> DupRedirNode:
    if redir.fd[0] == "var":
        redir.fd[1] = compile_command_argument(redir.fd[1])
    if redir.arg[0] == "var":
        redir.arg[1] = compile_command_argument(redir.arg[1])
    return redir


def compile_redirection_here(redir: HeredocRedirNode) -> HeredocRedirNode:
    if redir.fd[0] == "var":
        redir.fd[1] = compile_command_argument(redir.fd[1])
    redir.arg = compile_command_argument(redir.arg)
    return redir


def compile_redirection_single_arg(redir: SingleArgRedirNode) -> SingleArgRedirNode:
    if redir.fd[0] == "var":
        redir.fd[1] = compile_command_argument(redir.fd[1])
    return redir
