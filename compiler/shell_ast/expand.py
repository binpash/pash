import copy
import logging

from shasta.ast_node import *

################################################################################
# SAFE EXPANSION ANALYSIS
################################################################################

def log(msg: str):
    logging.info(f'Expansion: {msg}')

## This contains all necessary state of the expansion
class ExpansionState:
    variables: dict
    def __init__(self, variables: dict):
        self.variables = variables

    def __repr__(self):
        return f'ExpansionState: {self.variables}'

## This function checks if a word is safe to expand (i.e. if it will 
## not have unpleasant side-effects)
def safe_to_expand(arg_char: ArgChar):
    if isinstance(arg_char, VArgChar): # Variable
        return True
    return False

def guess_arg(arg):
    res = ""
    for arg_char in arg:
        if isinstance(arg_char, CArgChar) \
            or isinstance(arg_char, EArgChar):
            res += chr(arg_char.char)
        else:
            return None
    return res

def safe_arg(arg):
    return all([safe_arg_char(arg_char) for arg_char in arg])

def safe_args(args):
    return all([safe_arg(arg) for arg in args])

def safe_arg_char(arg_char: ArgChar):
    # character, escaped---noop, but safe
    if isinstance(arg_char, CArgChar) \
        or isinstance(arg_char, EArgChar):
        return True
    # tilde --- only reads system state, safe to do early assuming no writes to HOME prior
    elif isinstance(arg_char, TArgChar):
        return True # TODO 2020-11-24 MMG modified variable set? take in/output written vars...
    # arithmetic -- depends on what we have
    elif isinstance(arg_char, AArgChar):
        return safe_arith(arg_char.arg)
    # quoted -- safe if its contents are safe
    elif isinstance(arg_char, QArgChar):
        return safe_arg(arg_char.arg)
    # variables -- safe if the format is safe as are the remaining words
    elif isinstance(arg_char, VArgChar):
        return safe_var(fmt=arg_char.fmt, 
                        null=arg_char.null, 
                        var=arg_char.var, 
                        arg=arg_char.arg)
    # command substitution -- depends on the command
    elif isinstance(arg_char, BArgChar):
        return safe_command(arg_char.node)
    
    raise ValueError("bad object {}, expected one of CETAVQB".format(arg_char))

def safe_var(fmt, null, var, arg):
    if (fmt in ['Normal', 'Length']):
        return True
    elif (fmt in ['Minus', 'Plus', 'TrimR', 'TrimRMax', 'TrimL', 'TrimLMax']):
        return safe_arg(arg)
    elif (fmt in ['Question']):
        return False # TODO 2020-12-10 MMG unless we know `var` is set
    elif (fmt in ['Assign']):
        return False # TODO 2020-11-24 MMG unless we know `var` is set

    raise ValueError("bad parameter format {}".format(fmt))

def safe_arith(arg):
    # operations are safe
    # `+=` and `=` and family are UNSAFE
    # NONPOSIX: `++` and `--` are UNSAFE
    # `op="+=1"; $((x $op))` is UNSAFE

    # to determine safety, we:
    #   (a) check that every arg_char here is safe
    #   (b) pre-parse it symbolically well enough to ensure that no mutating operations occur
    expr = guess_arg(arg)

    if (arg is None):
        # TODO 2020-11-25 MMG symbolic pre-parse?
        return False
    elif ('=' in expr or '++' in expr or '--' in expr):
        # TODO 2020-11-25 MMG false negatives: ==, >=, <=
        return False
    else:
        # it's a concrete string that doesn't have mutation operations in it... go for it!
        return True

safe_cases = {
        "Pipe": (lambda:
                 lambda ast_node: safe_pipe(ast_node)),
        "Command": (lambda:
                    lambda ast_node: safe_simple(ast_node)),
        "And": (lambda:
                lambda ast_node: safe_and_or_semi(ast_node)),
        "Or": (lambda:
               lambda ast_node: safe_and_or_semi(ast_node)),
        "Semi": (lambda:
                 lambda ast_node: safe_and_or_semi(ast_node)),
        "Redir": (lambda:
                  lambda ast_node: safe_redir_subshell(ast_node)),
        "Subshell": (lambda:
                     lambda ast_node: safe_redir_subshell(ast_node)),
        "Background": (lambda:
                       lambda ast_node: safe_background(ast_node)),
        "Defun": (lambda:
                  lambda ast_node: safe_defun(ast_node)),
        "For": (lambda:
                  lambda ast_node: safe_for(ast_node)),
        "While": (lambda:
                  lambda ast_node: safe_while(ast_node)),
        "Case": (lambda:
                  lambda ast_node: safe_case(ast_node)),
        "If": (lambda:
                  lambda ast_node: safe_if(ast_node))
        }

def safe_command(command):
    # TODO 2020-11-24 MMG which commands are safe to run in advance?
    # TODO 2020-11-24 MMG how do we differentiate it being safe to do nested expansions?
    global safe_cases
    return ast_match(command, safe_cases)

def safe_pipe(node):
    return False

safe_commands = ["echo", ":"]

def safe_simple(node: CommandNode):
    # TODO 2020-11-25 check redirs, assignments

    if (len(node.arguments) <= 0):
        return True

    cmd = guess_arg(node.arguments[0])
    if (cmd is None or cmd not in safe_commands):
        return False
    else:
        return safe_args(node.arguments[1:])

def safe_and_or_semi(node):
    return False

def safe_redir_subshell(node):
    return False

def safe_background(node):
    return False

def safe_defun(node):
    return False

def safe_for(node):
    return False

def safe_while(node):
    return False

def safe_case(node):
    return False

def safe_if(node):
    return False

################################################################################
# EARLY EXPANSION
################################################################################

# General approach:
#
# - expand_* functions try to expand the AST
#   + words return a string when it works, raises when it doesn't
#     TODO MMG 2020-12-14 really should return (intermediate?) fields, not a single string
#   + commands just set the structural bits appropriately

# when early expansion detects an error
class EarlyError(RuntimeError):
    def __init__(self, arg):
        self.arg = arg

class StuckExpansion(RuntimeError):
    def __init__(self, reason, *info):
        self.reason = reason
        self.info = info

class ImpureExpansion(RuntimeError):
    def __init__(self, reason, *info):
        self.reason = reason
        self.info = info

class Unimplemented(RuntimeError):
    def __init__(self, msg, ast):
        self.msg = msg
        self.ast = ast

class InvalidVariable(RuntimeError):
    def __init__(self, var, reason):
        self.var = var
        self.reason = reason

def lookup_variable(var, exp_state):
    expanded_var = lookup_variable_inner(var, exp_state)
    
    return None, expanded_var

## Looksup a variable and flattens it if it is an array 
def lookup_variable_inner(varname, exp_state: ExpansionState):
    value = lookup_variable_inner_core(varname, exp_state)
    if value is not None and not isinstance(value, str):
        ## TODO: This is not handled at the moment (and it is unclear if it should be).
        ##
        ## This is only returned when we are in an array
        raise Unimplemented("Expanded value is not None or a string", (varname, value))
    return value

## Looks up the variable and if it is unset it raises an error
def lookup_variable_inner_core(varname, exp_state: ExpansionState):
    value = lookup_variable_inner_unsafe(varname, exp_state)
    if value is None and is_u_set(exp_state):
        raise StuckExpansion("-u is set and variable was unset", varname)
    return value


def lookup_variable_inner_unsafe(varname, exp_state: ExpansionState):
    ## TODO: Is it in there? If we have -u and it is in there.
    _type, value = exp_state.variables.get(varname, [None, None])
    return value

## This function checks if the -u flag is set
def is_u_set(exp_state: ExpansionState):
    value = lookup_variable_inner_unsafe('-', exp_state)
    # log(f'Previous set status is: {value}')
    return "u" in value


def invalidate_variable(var, reason, exp_state):
    exp_state.variables[var] = [None, InvalidVariable(var, reason)]
    return exp_state


def expand_args(args, exp_state, quoted = False):
    res = []
    for arg in args:
        new = expand_arg(arg, exp_state, quoted = quoted)

        # expanded! add the string in
        res.append(new)

    splitted_args = split_args(res, exp_state)

    return splitted_args

def split_args(args, exp_state):
    _, ifs = lookup_variable("IFS", exp_state)

    if ifs is None:
        ifs = "\n\t "

    ifs = [ord(c) for c in ifs]

    res = []
    for arg in args:
        cur = []

        for c in arg:
            if isinstance(c, CArgChar) and c.char in ifs:
                 # split!
                 if len(cur) > 0: # TODO(mmg): or if val isn't IFS whitespace
                     res.append(cur)
                 cur = []
            else:
                cur.append(c)

        if len(cur) > 0:
            res.append(cur)

    return res

def char_code(c) -> ArgChar:
    if c in "'\\\"()${}[]*?":
        return EArgChar(ord(c))
    else:
        return CArgChar(ord(c))

def expand_arg(arg_chars, exp_state, quoted = False):
    # log(f'expanding arg {arg_chars}")
    res = []
    for arg_char in arg_chars:
        new = expand_arg_char(arg_char, quoted, exp_state)

        if isinstance(new, str):
            res += [char_code(c) for c in list(new)]
        else:
            res.extend(new)

    return res

def expand_arg_char(arg_char: ArgChar, quoted, exp_state):
    if isinstance(arg_char, CArgChar):
        if arg_char.char in ['*', '?', '{', '}', '[', ']'] and not quoted:
            raise Unimplemented("globbing", arg_char)

        return [arg_char]
    elif isinstance(arg_char, EArgChar):
        ## 2021-09-15 MMG Just guessing here
        if arg_char.char in ['*', '?', '{', '}', '[', ']'] and not quoted:
            raise Unimplemented("globbing", arg_char)
        return [arg_char]
    elif isinstance(arg_char, TArgChar):
        val = arg_char.string
        if val is None or val == "" or val == "None":
            _type, val = lookup_variable("HOME", exp_state)

            if isinstance(val, InvalidVariable):
                raise StuckExpansion("HOME invalid for ~", arg_char, val)
            elif val is None:
                return "~"
            else:
                return val
        else:
            # TODO 2020-12-10 getpwnam
            raise Unimplemented("~ with prefix", arg_char)
    elif isinstance(arg_char, AArgChar):
        # TODO 2020-12-10 arithmetic parser and evaluator
        raise Unimplemented("arithmetic", arg_char)
    elif isinstance(arg_char, QArgChar):
        return [QArgChar(expand_arg(arg_char.arg, exp_state, quoted = True))]
        # return [['Q', expand_arg(arg_char.arg, exp_state, quoted = True)]]
    elif isinstance(arg_char, VArgChar):
        return expand_var(fmt=arg_char.fmt, 
                          null=arg_char.null, 
                          var=arg_char.var, 
                          arg=arg_char.arg, 
                          quoted=quoted, 
                          exp_state=exp_state)
    elif isinstance(arg_char, BArgChar):
        # TODO 2020-12-10 run commands?
        raise ImpureExpansion("command substitution", arg_char)
    else:
        raise Unimplemented("weird object", arg_char)

def expand_var(fmt, null, var, arg, quoted, exp_state):
    # TODO 2020-12-10 special variables

    _type, value = lookup_variable(var, exp_state)

    log(f'Var: {var} value: {value}')

    if isinstance(value, InvalidVariable):
        raise StuckExpansion("couldn't expand invalid variable", value)

    if fmt == 'Normal':
        if value is None:
            return ""
        else:
            return value
    elif fmt == 'Length':
        if value is None:
            return "0"
        else:
            return str(len(value))
    elif fmt == 'Minus':
        if value is None or (null and value == ""):
            return expand_arg(arg, exp_state, quoted = quoted)
        else:
            return value
    elif fmt == 'Assign':
        if value is None or (null and value == ""):
            raise ImpureExpansion("assignment format on unset/null variable", value, arg)
        else:
            return value
    elif fmt == 'Plus':
        if value is None or (null and value == ""):
            return ""
        else:
            return expand_arg(arg, exp_state, quoted = quoted)
    elif fmt == 'Question':
        if value is None or (null and value == ""):
            # TODO 2020-12-10 more context probably helpful here
            raise EarlyError(expand_arg(arg, exp_state, quoted = quoted))
        else:
            return value
    elif fmt in ['TrimR', 'TrimRMax', 'TrimL', 'TrimLMax']:
        # TODO need patterns
        raise Unimplemented("patterns", [fmt, null, var, arg])
    else:
        raise ValueError("bad parameter format {}".format(fmt))

expand_cases = {
        "Pipe": (lambda exp_state:
                 lambda ast_node: expand_pipe(ast_node, exp_state)),
        "Command": (lambda exp_state:
                    lambda ast_node: expand_simple(ast_node, exp_state)),
        "And": (lambda exp_state:
                lambda ast_node: expand_and_or_semi(ast_node, exp_state)),
        "Or": (lambda exp_state:
               lambda ast_node: expand_and_or_semi(ast_node, exp_state)),
        "Semi": (lambda exp_state:
                 lambda ast_node: expand_and_or_semi(ast_node, exp_state)),
        "Redir": (lambda exp_state:
                  lambda ast_node: expand_redir_subshell(ast_node, exp_state)),
        "Subshell": (lambda exp_state:
                     lambda ast_node: expand_redir_subshell(ast_node, exp_state)),
        "Background": (lambda exp_state:
                       lambda ast_node: expand_background(ast_node, exp_state)),
        "Defun": (lambda exp_state:
                  lambda ast_node: expand_defun(ast_node, exp_state)),
        "For": (lambda exp_state:
                  lambda ast_node: expand_for(ast_node, exp_state)),
        "While": (lambda exp_state:
                  lambda ast_node: expand_while(ast_node, exp_state)),
        "Case": (lambda exp_state:
                  lambda ast_node: expand_case(ast_node, exp_state)),
        "If": (lambda exp_state:
                  lambda ast_node: expand_if(ast_node, exp_state))
        }

def expand_command(command, exp_state: ExpansionState):
    # TODO 2020-11-24 MMG which commands are safe to run in advance?
    # TODO 2020-11-24 MMG how do we differentiate it being safe to do nested expansions?
    global expand_cases
    return ast_match(command, expand_cases, exp_state)

def expand_pipe(node, exp_state):
    for i, n in enumerate(node.items):
        # copy environment to simulate subshell (no outer effect)
        node.items[i] = expand_command(n, copy.deepcopy(exp_state))

    return node

def expand_simple(node, exp_state):
    # TODO 2020-11-25 MMG is this the order bash does?
    node.redir_list = expand_redir_list(node.redir_list, exp_state)

    if len(node.assignments) > 0:
        raise ImpureExpansion('assignment', node.assignments)
    
    node.arguments = expand_args(node.arguments, exp_state)

    return node

def expand_redir_list(redir_list, exp_state):
    for (i, r) in enumerate(redir_list):
        redir_list[i] = expand_redir(r, exp_state)

    return redir_list

def expand_redir(redirection: RedirectionNode, exp_state):
    file_arg = expand_arg(redirection.arg, exp_state)

    redirection.arg = file_arg
    return redirection

def expand_and_or_semi(node, exp_state):
    node.left_operand = expand_command(node.left_operand, exp_state)
    node.right_operand = expand_command(node.right_operand, exp_state)

    return node

def expand_redir_subshell(node, exp_state):
    # copy environment to simulate subshell (no outer effect)
    node.node = expand_command(node.node, copy.deepcopy(exp_state))

    return node

def expand_background(node, exp_state):
    # copy environment to simulate subshell (no outer effect)
    node.node = expand_command(node.node, copy.deepcopy(exp_state))

    return node

def expand_defun(node, exp_state):
    # TODO 2020-11-24 MMG invalidate postional args
    node.body = expand_command(node.body, copy.deepcopy(exp_state))

    return node

def expand_for(node, exp_state):
    node.argument = expand_arg(node.argument, exp_state)

    # TODO 2020-11-24 if node.argument is fully expanded, we can just unroll the loop
    exp_state = invalidate_variable(node.variable, "variable of for loop", exp_state)
    node.body = expand_command(node.body, exp_state)

    return node

def expand_while(node, exp_state):
    node.test = expand_command(node.test, exp_state)
    node.body = expand_command(node.body, exp_state)

    return node

def expand_case(node, exp_state):
    # TODO 2020-11-24 preprocess scrutinee, each pattern, each case

    raise Unimplemented("case statements", node)

def expand_if(node, exp_state):
    node.cond = expand_command(node.cond, exp_state)
    node.then_b = expand_command(node.then_b, exp_state)
    node.else_b = expand_command(node.else_b, exp_state)

    return node
