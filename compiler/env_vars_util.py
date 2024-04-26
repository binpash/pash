import shlex
from datetime import datetime

from util import log, print_time_delta

def read_vars_file(var_file_path):
    log("Reading variables from:", var_file_path)

    if(not var_file_path is None):
        vars_dict = {}
        # with open(var_file_path) as f:
        #     lines = [line.rstrip() for line in f.readlines()]

        with open(var_file_path, 'rb') as f:
            variable_reading_start_time = datetime.now()
            data = f.read().decode('utf-8', errors='replace')
            variable_reading_end_time = datetime.now()
            print_time_delta("Variable Reading", variable_reading_start_time, variable_reading_end_time)

            variable_tokenizing_start_time = datetime.now()
            ## TODO: Can we replace this tokenizing process with our own code? This is very slow :'(
            ##       It takes about 15ms on deathstar.
            tokens = shlex.split(data)
            variable_tokenizing_end_time = datetime.now()
            print_time_delta("Variable Tokenizing", variable_tokenizing_start_time, variable_tokenizing_end_time)
            # log("Tokens:", tokens)

        # MMG 2021-03-09 definitively breaking on newlines (e.g., IFS) and function outputs (i.e., `declare -f`)
        # KK  2021-10-26 no longer breaking on newlines (probably)

        ## At the start of each iteration token_i should point to a 'declare'
        token_i = 0
        while token_i < len(tokens):
            # FIXME is this assignment needed?
            export_or_typeset = tokens[token_i]

            ## Array variables require special parsing treatment
            if (export_or_typeset == "declare" and is_array_variable(tokens[token_i+1])):
                var_name, var_type, var_value, new_token_i = parse_array_variable(tokens, token_i)
                vars_dict[var_name] = (var_type, var_value)
                token_i = new_token_i
                continue

            new_token_i = find_next_delimiter(tokens, token_i)
            rest = " ".join(tokens[(token_i+1):new_token_i])
            token_i = new_token_i

            space_index = rest.find(' ')
            eq_index = rest.find('=')
            var_type = None

            ## Declared but unset?
            if eq_index == -1:
                if space_index != -1:
                    var_name = rest[(space_index+1):]
                    var_type = rest[:space_index]
                else:
                    var_name = rest
                var_value = ""
            ## Set, with type
            elif(space_index < eq_index and not space_index == -1):
                var_type = rest[:space_index]

                if var_type == "--":
                    var_type = None
                
                var_name = rest[(space_index+1):eq_index]
                var_value = rest[(eq_index+1):]
            ## Set, without type
            else:
                var_name = rest[:eq_index]
                var_value = rest[(eq_index+1):]

            ## Strip quotes
            if var_value is not None and len(var_value) >= 2 and \
               var_value[0] == "\"" and var_value[-1] == "\"":
                var_value = var_value[1:-1]                
                
            vars_dict[var_name] = (var_type, var_value)

        final_vars_dict = set_special_parameters(vars_dict)
        return final_vars_dict


## This sets the values of the special shell parameters correctly
##
## TODO KK PR#246 Do we need to split using IFS or is it always spaces?
##
## TODO MMG this isn't quite adequate: if pash_input_args contains
##      spaces, we'll miscount. KK and I wrote a test
##      evaluation/tests/interface_tests that's disabled as of PR#246.
##
##      the right solution here is:
##
##         - positional arguments get their own field in the
##           exp_state---they're not store with ordinary shell
##           variables
##
##         - we save those separately, probably in a separate file
##
##           ```
##           echo pash_argc=$# >pash_positional_args
##           for i in $(seq 0 $#)
##           do
##             echo "pash_arg$i=\"$i\"" >pash_positional_args
##           done
##           ```
##
##         - we load these separately. pretty annoying; here's a sketch
##
##           ```
##           cmd="set --"
##           for i in $(seq 0 $pash_argc)
##           do
##             cmd="$cmd \"\$pash_arg$i\""
##           done
##           eval "$cmd"
def set_special_parameters(variables: dict):
    new_vars = variables.copy()

    ia_t, input_args = get_var(variables, 'pash_input_args')
    es_t, exit_status = get_var(variables, 'pash_previous_exit_status')
    ss_t, set_status = get_var(variables, 'pash_previous_set_status')
    sn_t, shell_name = get_var(variables, 'pash_shell_name')

    ## TODO: Set the types of variables correctly
    new_vars['@'] = ia_t, " ".join(input_args)
    new_vars['?'] = es_t, exit_status
    new_vars['-'] = ss_t, set_status
    new_vars['0'] = sn_t, shell_name
    new_vars['#'] = ia_t, str(len(input_args))

    for i, arg in enumerate(input_args):
        index = i + 1
        new_vars[str(index)] = input_args[i]

    return new_vars

def get_var(variables: dict, varname: str):
    type, value = variables.get(varname, [None, None])
    return type, value

def is_array_variable(token):
    return ('a' in token)

## Based on the following:
## https://www.gnu.org/software/bash/manual/html_node/ANSI_002dC-Quoting.html#ANSI_002dC-Quoting
def ansi_c_expand(string):
    return bytes(string, "utf-8").decode("unicode_escape")

## This finds the end of this variable/function
def find_next_delimiter(tokens, i):
    if (tokens[i] == "declare"):
        return i + 3
    else:
        ## TODO: When is this case actually useful?
        j = i + 1
        while j < len(tokens) and (tokens[j] != "declare"):
            j += 1
        return j

def parse_array_variable(tokens, i):
    ## The `declare` keyword
    _declare = tokens[i]
    ## The type
    declare_type = tokens[i+1]
    assert(is_array_variable(declare_type))

    ## The variable name and first argument
    ## TODO: Test with empty array and single value array
    name_and_start=tokens[i+2]
    first_equal_index = name_and_start.find('=')

    ## If it doesn't contain any = then it is empty
    if first_equal_index == -1:
        ## Then the name is the whole token,
        ##  the type is None (TODO)
        ##  and the value is empty
        return name_and_start, None, "", i+3

    var_name = name_and_start[:first_equal_index]
    array_start = name_and_start[first_equal_index+1:]

    var_values = []
    if array_start == "()":
        next_i = i+3
    else:
        ## Remove the opening parenthesis
        array_item = array_start[1:]

        ## Set the index that points to array items
        curr_i = i+2

        done = False
        while not done:
            ## TODO: Is this check adequate? Or could it miss the end 
            ##       (or be misleaded into an earlier end by the item value?)
            if array_item.endswith(")"):
                done = True
                array_item = array_item[:-1]

            first_equal_index = array_item.find('=')
            ## Find the index and value of the array item
            item_index_raw = array_item[:first_equal_index]
            item_value = array_item[first_equal_index+1:]

            ## Sometimes the value starts with a dollar mark, see Bash ANSI-C quoting:
            ## https://www.gnu.org/software/bash/manual/html_node/ANSI_002dC-Quoting.html#ANSI_002dC-Quoting
            if item_value.startswith("$"):
                ## TODO: Figure out if this is adequate
                item_value = ansi_c_expand(item_value[1:])

            item_index = int(item_index_raw[1:-1])
            
            ## Add None values if the index is larger than the next item (see Bash sparse arrays)
            ## TODO: Keep bash array values as maps to avoid sparse costs 
            var_values += [None] * (item_index - len(var_values))
            ## Set the next item
            var_values.append(item_value)

            

            ## Get next array_item
            curr_i += 1
            array_item = tokens[curr_i]
        
        next_i = curr_i

    ## TODO: Michael?
    var_type = None

    return var_name, var_type, var_values, next_i
