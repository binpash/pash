from shasta import ast_node
import pexpect
import config
from util import *

def expand_using_bash(ast: ast_node) -> ast_node:
    pass

def wait_bash_mirror(bash_mirror):
    r = bash_mirror.expect(r'EXPECT\$ ')
    assert(r == 0)
    output = bash_mirror.before

    ## I am not sure why, but \r s are added before \n s
    output = output.replace('\r\n', '\n')

    log("Before the prompt!")
    log(output)
    return output


def query_expand_variable_bash_mirror(variable):
    global bash_mirror

    command = f'if [ -z ${{{variable}+foo}} ]; then echo -n "PASH_VAR_UNSET"; else echo -n "${variable}"; fi'
    data = sync_run_line_command_mirror(command)

    if data == "PASH_VAR_UNSET":
        return None
    else:
        ## This is here because we haven't specified utf encoding when spawning bash mirror
        # return data.decode('ascii')
        return data

def query_expand_bash_mirror(string):
    global bash_mirror

    command = f'echo -n "{string}"'
    return sync_run_line_command_mirror(command)

def sync_run_line_command_mirror(command):
    bash_command = f'{command}'
    log("Executing bash command in mirror:", bash_command)

    bash_mirror.sendline(bash_command)

    data = wait_bash_mirror(bash_mirror)
    log("mirror done!")

    return data


def update_bash_mirror_vars(var_file_path):
    global bash_mirror

    assert(var_file_path != ""  and not var_file_path is None)

    bash_mirror.sendline(f'PS1="EXPECT\$ "')
    wait_bash_mirror(bash_mirror)
    log("PS1 set!")

    ## TODO: There is unnecessary write/read to this var file now.
    bash_mirror.sendline(f'source {var_file_path}')
    log("sent source to mirror")
    wait_bash_mirror(bash_mirror)
    log("mirror done!")

def add_to_variable_cache(variable_name, value):
    global variable_cache
    variable_cache[variable_name] = value

def get_from_variable_cache(variable_name):
    global variable_cache
    try:
        return variable_cache[variable_name]
    except:
        return None

def reset_variable_cache():
    global variable_cache

    variable_cache = {}


## This finds the end of this variable/function
def find_next_delimiter(tokens, i):
    if (tokens[i] == "declare"):
        return i + 3
    else:
        j = i + 1
        while j < len(tokens) and (tokens[j] != "declare"):
            j += 1
        return j
