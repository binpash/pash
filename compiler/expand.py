from typing import NoReturn
from definitions.ast_node import *
from definitions.ast_node_c import *

import ast_to_ir
import config

################################################################################
# SAFE EXPANSION ANALYSIS
################################################################################

## This function checks if a word is safe to expand (i.e. if it will 
## not have unpleasant side-effects)
def safe_to_expand(arg_char):
    key, val = get_kv(arg_char)
    if (key in ['V']): # Variable
        return True
    return False

def guess_arg(arg):
    res = ""
    for arg_char in arg:
        key, val = get_kv(arg_char)

        if (key in ['C', 'E']):
            res += chr(val)
        else:
            return None
    return res

def safe_arg(arg):
    return all([safe_arg_char(arg_char) for arg_char in arg])

def safe_args(args):
    return all([safe_arg(arg) for arg in args])

def safe_arg_char(arg_char):
    key, val = get_kv(arg_char)
    # character, escaped---noop, but safe
    if (key in ['C', 'E']): 
        return True
    # tilde --- only reads system state, safe to do early assuming no writes to HOME prior
    elif (key == 'T'):
        return True # TODO 2020-11-24 MMG modified variable set? take in/output written vars...
    # arithmetic -- depends on what we have
    elif (key == 'A'):
        return safe_arith(val)
    # quoted -- safe if its contents are safe
    elif (key == 'Q'):
        return safe_arg(val)
    # variables -- safe if the format is safe as are the remaining words
    elif (key == 'V'):
        return safe_var(*val)
    # command substitution -- depends on the command
    elif (key == 'B'):
        return safe_command(val)
    
    raise ValueError("bad key {}, expected one of CETAVQB".format(key))

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
    return ast_to_ir.ast_match(command, safe_cases)

def safe_pipe(node):
    return False

safe_commands = ["echo", ":"]

def safe_simple(node):
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
#   + words return a string when it works, None when it doesn't
#   + commands just set the structural bits appropriately

# when early expansion detects an error
class EarlyError(RuntimeError):
    def __init__(self, arg):
        self.arg = arg

def expand_args(args, quoted, config):
    res = []
    for arg in args:
        new = expand_arg(arg, quoted, config)

        if new is None:
            # expansion failed; just keep rolling
            res.append(arg)
        else:
            # expanded! add the string in
            res.append(new)

    return res

def expand_arg(arg_chars, quoted, config):
    res = []
    for arg_char in arg_chars:
        new = expand_arg_char(arg_char, quoted, config)

        if isinstance(new, str):
            res += [["E", ord(c)] for c in list(new)]
        elif new is None:
            res.append(arg_char)
        else:
            res.extend(new)

    return res

def expand_arg_char(arg_char, quoted, config):
    key, val = get_kv(arg_char)

    if key in ['C', 'E']:
        return
    elif key == 'T':
        if val == "":
            if "HOME" in config['shell_variables']:
                return config['shell_variables']["HOME"]
            else:
                return "~"
        else:
            # TODO 2020-12-10 getpwnam
            return
    elif key == 'A':
        # TODO 2020-12-10 arithmetic parser and evaluator
        return
    elif key == 'Q':
        return expand_arg(val, True, config)
    elif key == 'V':
        fmt, null, var, arg = val
        return expand_var(fmt, null, var, arg, quoted, config)
    elif key == 'B':
        # TODO 2020-12-10 run commands?
        return

def expand_var(fmt, null, var, arg, quoted, config):
    # TODO 2020-12-10 special variables

    if var in config['shell_variables']:
        type, value = config['shell_variables'][var]
    else:
        type = None
        value = None

    if fmt == 'Normal':
        return value
    elif fmt == 'Length':
        if value is None:
            return "0"
        else:
            return str(len(value))
    elif fmt == 'Minus':
        if value is None or (null and value == ""):
            return expand_arg(arg, quoted, config)
        else:
            return value
    elif fmt == 'Assign':
        if value is None or (null and value == ""):
            new = expand_arg(arg, quoted, config)
            if new is not None:
                config['shell_variables'][var] = (None, new)
                return new
            else:
                return
        else:
            return value
    elif fmt == 'Plus':
        if value is None or (null and value == ""):
            return ""
        else:
            return expand_arg(arg, quoted, config)
    elif fmt == 'Question':
        if value is None or (null and value == ""):
            # TODO 2020-12-10 more context probably helpful here
            raise EarlyError(expand_arg(arg, quoted, config))
        else:
            return value
    elif fmt == 'TrimR':
        # TODO need patterns
        return
    elif fmt == 'TrimRMax':
        return
    elif fmt == 'TrimL':
        return
    elif fmt == 'TrimLMax':
        return
    else:
        raise ValueError("bad parameter format {}".format(fmt))

expand_cases = {
        "Pipe": (lambda config:
                 lambda ast_node: expand_pipe(ast_node, config)),
        "Command": (lambda config:
                    lambda ast_node: expand_simple(ast_node, config)),
        "And": (lambda config:
                lambda ast_node: expand_and_or_semi(ast_node, config)),
        "Or": (lambda config:
               lambda ast_node: expand_and_or_semi(ast_node, config)),
        "Semi": (lambda config:
                 lambda ast_node: expand_and_or_semi(ast_node, config)),
        "Redir": (lambda config:
                  lambda ast_node: expand_redir_subshell(ast_node, config)),
        "Subshell": (lambda config:
                     lambda ast_node: expand_redir_subshell(ast_node, config)),
        "Background": (lambda config:
                       lambda ast_node: expand_background(ast_node, config)),
        "Defun": (lambda config:
                  lambda ast_node: expand_defun(ast_node, config)),
        "For": (lambda config:
                  lambda ast_node: expand_for(ast_node, config)),
        "While": (lambda config:
                  lambda ast_node: expand_while(ast_node, config)),
        "Case": (lambda config:
                  lambda ast_node: expand_case(ast_node, config)),
        "If": (lambda config:
                  lambda ast_node: expand_if(ast_node, config))
        }

def expand_command(command, config):
    # TODO 2020-11-24 MMG which commands are safe to run in advance?
    # TODO 2020-11-24 MMG how do we differentiate it being safe to do nested expansions?
    global expand_cases
    return ast_to_ir.ast_match(command, expand_cases, config)

def expand_pipe(node, config):
    return False

def expand_simple(node, config):
    # TODO 2020-11-25 check redirs, assignments

    if (len(node.arguments) == 0):
        return node

    node.arguments = expand_args(node.arguments, False, config)

    return node

def expand_and_or_semi(node, config):
    return node

def expand_redir_subshell(node, config):
    return node

def expand_background(node, config):
    return node

def expand_defun(node, config):
    return node

def expand_for(node, config):
    return node

def expand_while(node, config):
    return node

def expand_case(node, config):
    return node

def expand_if(node, config):
    return node