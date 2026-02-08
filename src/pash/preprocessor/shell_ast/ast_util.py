"""
AST utility functions for the PaSh preprocessor.
"""

from env_var_names import loop_iters_var, loop_iter_var
from shasta.ast_node import AstNode, ArgChar, CArgChar
from util import make_kv, unzip


class PreprocessedAST:
    """Class used by the preprocessor in ast_to_ir."""

    def __init__(
        self, ast, replace_whole, non_maximal, something_replaced=True, last_ast=False
    ):
        assert isinstance(ast, AstNode)
        self.ast = ast
        self.replace_whole = replace_whole
        self.non_maximal = non_maximal
        self.something_replaced = something_replaced
        self.last_ast = last_ast

    def should_replace_whole_ast(self):
        return self.replace_whole

    def is_non_maximal(self):
        return self.non_maximal

    def will_anything_be_replaced(self):
        return self.something_replaced

    def is_last_ast(self):
        return self.last_ast


class UnparsedScript:
    """
    Represents text that was not modified at all by preprocessing,
    and therefore does not need to be unparsed.
    """

    def __init__(self, text):
        self.text = text


# === Pattern matching for the AST ===


def format_args(args):
    formatted_args = [format_arg_chars(arg_chars) for arg_chars in args]
    return formatted_args


def format_arg_chars(arg_chars):
    chars = [format_arg_char(arg_char) for arg_char in arg_chars]
    return "".join(chars)


def format_arg_char(arg_char: ArgChar) -> str:
    return arg_char.format()


def string_to_carg_char_list(string: str) -> "list[CArgChar]":
    ret = [CArgChar(ord(char)) for char in string]
    return ret


def string_to_arguments(string):
    return [string_to_argument(word) for word in string.split(" ")]


def string_to_argument(string):
    ret = [char_to_arg_char(char) for char in string]
    return ret


def concat_arguments(arg1, arg2):
    """Arguments are simply `arg_char list` and therefore can just be concatenated."""
    return arg1 + arg2


def char_to_arg_char(char):
    return ["C", ord(char)]


def escaped_char(char):
    return ["E", ord(char)]


def standard_var_ast(string):
    return make_kv("V", ["Normal", False, string, []])


def make_arith(arg):
    return make_kv("A", arg)


def make_quoted_variable(string):
    return make_kv("Q", [standard_var_ast(string)])


def quote_arg(arg):
    return make_kv("Q", arg)


def redir_append_stderr_to_string_file(string):
    return make_kv("File", ["Append", 2, string_to_argument(string)])


def redir_stdout_to_file(arg):
    return make_kv("File", ["To", 1, arg])


def redir_file_to_stdin(arg):
    return make_kv("File", ["From", 0, arg])


def make_background(body, redirections=None):
    redirections = [] if redirections is None else redirections
    lineno = 0
    node = make_kv("Background", [lineno, body, redirections])
    return node


def make_backquote(node):
    node = make_kv("B", node)
    return node


def make_subshell(body, redirections=None):
    redirections = [] if redirections is None else redirections
    lineno = 0
    node = make_kv("Subshell", [lineno, body, redirections])
    return node


def make_command(arguments, redirections=None, assignments=None):
    redirections = [] if redirections is None else redirections
    assignments = [] if assignments is None else assignments
    lineno = 0
    node = make_kv("Command", [lineno, assignments, arguments, redirections])
    return node


def make_nop():
    return make_command([string_to_argument(":")])


def make_assignment(var, value):
    lineno = 0
    assignment = (var, value)
    assignments = [assignment]
    node = make_kv("Command", [lineno, assignments, [], []])
    return node


def make_semi_sequence(asts):
    if len(asts) == 0:
        return make_nop()

    if len(asts) == 1:
        return asts[0]
    else:
        acc = asts[-1]
        # Remove the last ast
        iter_asts = asts[:-1]
        for ast in iter_asts[::-1]:
            acc = make_kv("Semi", [ast, acc])
        return acc


def make_defun(name, body):
    lineno = 0
    node = make_kv("Defun", [lineno, name, body])
    return node


# === Make some nodes ===


def make_export_var_constant_string(var_name: str, value: str):
    node = make_export_var(var_name, string_to_argument(value))
    return node


def make_export_var(var_name: str, arg_char_list):
    """An argument is an arg_char_list."""
    arg1 = string_to_argument(f"{var_name}=")
    arguments = [string_to_argument("export"), concat_arguments(arg1, arg_char_list)]
    node = make_command(arguments)
    return node


def export_pash_loop_iters_for_current_context(all_loop_ids: "list[int]"):
    if len(all_loop_ids) > 0:
        iter_var_names = [loop_iter_var(loop_id) for loop_id in all_loop_ids]
        iter_vars = [
            standard_var_ast(iter_var_name) for iter_var_name in iter_var_names
        ]
        concatted_vars = [iter_vars[0]]
        for iter_var in iter_vars[1:]:
            concatted_vars.append(char_to_arg_char("-"))
            concatted_vars.append(iter_var)
        quoted_vars = [quote_arg(concatted_vars)]
    else:
        quoted_vars = []

    # export pash_loop_iters="$pash_loop_XXX_iter $pash_loop_YYY_iter ..."
    save_loop_iters_node = make_export_var(loop_iters_var(), quoted_vars)

    return save_loop_iters_node


def make_unset_var(var_name: str):
    arguments = [string_to_argument("unset"), string_to_argument(var_name)]
    node = make_command(arguments)
    return node


def make_loop_list_assignment(loop_list_args):
    """Create HS_LOOP_LIST=<list> assignment from for-loop arguments.

    Matches the implementation from fae47999 commit of spec_future branch.

    Args:
        loop_list_args: list[list[ArgChar]] - the for-loop's argument field

    Returns:
        AST node for HS_LOOP_LIST=<list> assignment
    """
    from shasta.ast_node import CArgChar, QArgChar, to_ast_node
    import copy

    # Create initial assignment node with placeholder value '0'
    list_eval_node = to_ast_node(make_assignment('HS_LOOP_LIST', string_to_argument('0')))

    # Concatenate all for-loop arguments with spaces (matching fae47999 implementation)
    list_arguments = copy.deepcopy(loop_list_args[0])
    for a in loop_list_args[1:]:
        list_arguments.append(CArgChar(ord(' ')))
        list_arguments.extend(copy.deepcopy(a))

    # Set the assignment value to the quoted list (QArgChar wrapper)
    list_eval_node.assignments[0].val = [QArgChar(list_arguments)]

    return list_eval_node


def make_increment_var(var_name: str):
    arg = string_to_argument(f"{var_name}+1")
    arith_expr = make_arith(arg)
    assignments = [[var_name, [arith_expr]]]
    node = make_command([], assignments=assignments)
    return node
