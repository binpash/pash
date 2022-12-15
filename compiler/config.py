import json
import os
import subprocess
import math
import shlex

from datetime import datetime

from ir_utils import *
from util import *

## Global
__version__ = "0.10" # FIXME add libdash version
GIT_TOP_CMD = [ 'git', 'rev-parse', '--show-toplevel', '--show-superproject-working-tree']
if 'PASH_TOP' in os.environ:
    PASH_TOP = os.environ['PASH_TOP']
else:
    PASH_TOP = subprocess.run(GIT_TOP_CMD, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True).stdout.rstrip()

PYTHON_VERSION = "python3"
PLANNER_EXECUTABLE = os.path.join(PASH_TOP, "compiler/pash_runtime.py")
RUNTIME_EXECUTABLE = os.path.join(PASH_TOP, "compiler/pash_runtime.sh")
SAVE_ARGS_EXECUTABLE = os.path.join(PASH_TOP, "runtime/save_args.sh")

## Ensure that PASH_TMP_PREFIX is set by pa.sh
assert(not os.getenv('PASH_TMP_PREFIX') is None)
PASH_TMP_PREFIX = os.getenv('PASH_TMP_PREFIX')

##
## Global configuration used by all pash components
##
LOGGING_PREFIX = ""
OUTPUT_TIME = False
DEBUG_LEVEL = 0
LOG_FILE = ""


HDFS_PREFIX = "$HDFS_DATANODE_DIR/"


config = {}
annotations = []
pash_args = None


## This function sets the global configuration
##
## TODO: Actually move everything outside of pash_args to configuration.
def set_config_globals_from_pash_args(given_pash_args):
    global pash_args, OUTPUT_TIME, DEBUG_LEVEL, LOG_FILE
    pash_args = given_pash_args
    OUTPUT_TIME = pash_args.output_time
    DEBUG_LEVEL = pash_args.debug
    LOG_FILE = pash_args.log_file


## Increase the recursion limit (it seems that the parser/unparser needs it for bigger graphs)
sys.setrecursionlimit(10000)

def load_config(config_file_path=""):
    global config
    pash_config = {}
    CONFIG_KEY = 'distr_planner'

    if(config_file_path == ""):
      config_file_path = '{}/compiler/config.json'.format(PASH_TOP)
    with open(config_file_path) as config_file:
        pash_config = json.load(config_file)

    if not pash_config:
        raise Exception('No valid configuration could be loaded from {}'.format(config_file_path))

    if CONFIG_KEY not in pash_config:
        raise Exception('Missing `{}` config in {}'.format(CONFIG_KEY, config_file_path))

    config = pash_config

def getWidth():
    cpus = os.cpu_count()
    return math.floor(cpus / 8) if cpus >= 16 else 2

def add_general_config_arguments(parser):
    ## TODO: Delete that at some point, or make it have a different use (e.g., outputting time even without -d 1).
    parser.add_argument("-t", "--output_time", #FIXME: --time
                        help="(obsolete, time is always logged now) output the time it took for every step",
                        action="store_true")
    parser.add_argument("-d", "--debug",
                        type=int,
                        help="configure debug level; defaults to 0",
                        default=0)
    parser.add_argument("--log_file",
                        help="configure where to write the log; defaults to stderr.",
                        default="")

## These are arguments that are common to pash.py and pash_runtime.py
def add_common_arguments(parser):
    add_general_config_arguments(parser)

    parser.add_argument("-w", "--width",
                        type=int,
                        default=getWidth(),
                        help="set data-parallelism factor")
    parser.add_argument("--no_optimize",
                        help="not apply transformations over the DFG",
                        action="store_true")
    parser.add_argument("--dry_run_compiler",
                        help="not execute the compiled script, even if the compiler succeeded",
                        action="store_true")
    parser.add_argument("--assert_compiler_success",
                        help="assert that the compiler succeeded (used to make tests more robust)",
                        action="store_true")
    parser.add_argument("--avoid_pash_runtime_completion",
                        help="avoid the pash_runtime execution completion (only relevant when --debug > 0)",
                        action="store_true")
    parser.add_argument("--profile_driven",
                        help="(experimental) use profiling information when optimizing",
                        action="store_true")
    parser.add_argument("-p", "--output_optimized", # FIXME: --print
                        help="output the parallel shell script for inspection",
                        action="store_true")
    parser.add_argument("--graphviz",
                        help="generates graphical representations of the dataflow graphs. The option argument corresponds to the format. PaSh stores them in a timestamped directory in the argument of --graphviz_dir",
                        choices=["no", "dot", "svg", "pdf", "png"],
                        default="no")
    ## TODO: To discuss: Do we maybe want to have graphviz to always be included 
    ##       in the temp directory (under a graphviz subdirectory) instead of in its own?
    ##   kk: I think that ideally we want a log-directory where we can put logs, graphviz, 
    ##       and other observability and monitoring info (instead of putting them in the temp).
    parser.add_argument("--graphviz_dir",
                        help="the directory in which to store graphical representations",
                        default="/tmp")
    parser.add_argument("--no_eager",
                        help="(experimental) disable eager nodes before merging nodes",
                        action="store_true")
    parser.add_argument("--no_cat_split_vanish",
                        help="(experimental) disable the optimization that removes cat with N inputs that is followed by a split with N inputs",
                        action="store_true")
    parser.add_argument("--no_daemon",
                        help="Run the compiler everytime we need a compilation instead of using the daemon",
                        action="store_true",
                        default=False)
    parser.add_argument("--parallel_pipelines",
                        help="Run multiple pipelines in parallel if they are safe to run",
                        action="store_true",
                        default=False)
    parser.add_argument("--r_split",
                        help="(experimental) use round robin split, merge, wrap, and unwrap",
                        action="store_true")
    parser.add_argument("--r_split_batch_size",
                        type=int,
                        help="(experimental) configure the batch size of r_split (default: 1MB)",
                        default=1000000)
    parser.add_argument("--dgsh_tee",
                        help="(experimental) use dgsh-tee instead of eager",
                        action="store_true")
    parser.add_argument("--speculation",
                        help="(experimental) run the original script during compilation; if compilation succeeds, abort the original and run only the parallel (quick_abort) (Default: no_spec)",
                        choices=['no_spec', 'quick_abort'],
                        default='no_spec')
    parser.add_argument("--termination",
                        help="(experimental) determine the termination behavior of the DFG. Defaults to cleanup after the last process dies, but can drain all streams until depletion",
                        choices=['clean_up_graph', 'drain_stream'],
                        default="clean_up_graph")
    parser.add_argument("--daemon_communicates_through_unix_pipes",
                        help="(experimental) the daemon communicates through unix pipes instead of sockets",
                        action="store_true")
    parser.add_argument("--distributed_exec",
                        help="(experimental) execute the script in a distributed environment. Remote machines should be configured and ready",
                        action="store_true",
                        default=False)
    parser.add_argument("--config_path",
                        help="determines the config file path. By default it is 'PASH_TOP/compiler/config.yaml'.",
                        default="")
    parser.add_argument("--version",
            action='version',
            version='%(prog)s {version}'.format(version=__version__))
    return

def pass_common_arguments(pash_arguments):
    arguments = []
    if (pash_arguments.no_optimize):
        arguments.append(string_to_argument("--no_optimize"))
    if (pash_arguments.dry_run_compiler):
        arguments.append(string_to_argument("--dry_run_compiler"))
    if (pash_arguments.assert_compiler_success):
        arguments.append(string_to_argument("--assert_compiler_success"))
    if (pash_arguments.avoid_pash_runtime_completion):
        arguments.append(string_to_argument("--avoid_pash_runtime_completion"))
    if (pash_arguments.profile_driven):
        arguments.append(string_to_argument("--profile_driven"))
    if (pash_arguments.output_time):
        arguments.append(string_to_argument("--output_time"))
    if (pash_arguments.output_optimized):
        arguments.append(string_to_argument("--output_optimized"))
    arguments.append(string_to_argument("--graphviz"))
    arguments.append(string_to_argument(pash_arguments.graphviz))
    arguments.append(string_to_argument("--graphviz_dir"))
    arguments.append(string_to_argument(pash_arguments.graphviz_dir))
    if(not pash_arguments.log_file == ""):
        arguments.append(string_to_argument("--log_file"))
        arguments.append(string_to_argument(pash_arguments.log_file))
    if (pash_arguments.no_eager):
        arguments.append(string_to_argument("--no_eager"))
    if (pash_arguments.r_split):
        arguments.append(string_to_argument("--r_split"))
    if (pash_arguments.dgsh_tee):
        arguments.append(string_to_argument("--dgsh_tee"))
    if (pash_arguments.no_daemon):
        arguments.append(string_to_argument("--no_daemon"))
    if (pash_arguments.distributed_exec):
        arguments.append(string_to_argument("--distributed_exec"))
    if (pash_arguments.parallel_pipelines):
        arguments.append(string_to_argument("--parallel_pipelines"))
    if (pash_arguments.daemon_communicates_through_unix_pipes):
        arguments.append(string_to_argument("--daemon_communicates_through_unix_pipes"))
    arguments.append(string_to_argument("--r_split_batch_size"))
    arguments.append(string_to_argument(str(pash_arguments.r_split_batch_size)))
    if (pash_arguments.no_cat_split_vanish):
        arguments.append(string_to_argument("--no_cat_split_vanish"))
    arguments.append(string_to_argument("--debug"))
    arguments.append(string_to_argument(str(pash_arguments.debug)))
    arguments.append(string_to_argument("--termination"))
    arguments.append(string_to_argument(pash_arguments.termination))
    arguments.append(string_to_argument("--speculation"))
    arguments.append(string_to_argument(pash_arguments.speculation))
    arguments.append(string_to_argument("--width"))
    arguments.append(string_to_argument(str(pash_arguments.width)))
    if(not pash_arguments.config_path == ""):
        arguments.append(string_to_argument("--config_path"))
        arguments.append(string_to_argument(pash_arguments.config_path))
    return arguments

def init_log_file():
    global LOG_FILE
    if(not LOG_FILE == ""):
        with open(LOG_FILE, "w") as f:
            pass


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

##
## Read a shell variables file
##

def read_vars_file(var_file_path):
    global config

    log("Reading variables from:", var_file_path)


    config['shell_variables'] = None
    config['shell_variables_file_path'] = var_file_path
    if(not var_file_path is None):
        vars_dict = {}
        # with open(var_file_path) as f:
        #     lines = [line.rstrip() for line in f.readlines()]

        with open(var_file_path) as f:
            variable_reading_start_time = datetime.now()
            data = f.read()
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

        config['shell_variables'] = vars_dict
