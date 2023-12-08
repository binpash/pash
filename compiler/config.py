import json
import logging
import os
import subprocess
import math

from util import *

## Global
__version__ = "0.12.2"  # FIXME add libdash version
GIT_TOP_CMD = [
    "git",
    "rev-parse",
    "--show-toplevel",
    "--show-superproject-working-tree",
]
if "PASH_TOP" in os.environ:
    PASH_TOP = os.environ["PASH_TOP"]
else:
    PASH_TOP = subprocess.run(
        GIT_TOP_CMD,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
    ).stdout.rstrip()

PYTHON_VERSION = "python3"
PLANNER_EXECUTABLE = os.path.join(PASH_TOP, "compiler/pash_compiler.py")
RUNTIME_EXECUTABLE = os.path.join(PASH_TOP, "compiler/pash_runtime.sh")
SAVE_ARGS_EXECUTABLE = os.path.join(PASH_TOP, "runtime/save_args.sh")
SAVE_SHELL_STATE_EXECUTABLE = os.path.join(
    PASH_TOP, "compiler/orchestrator_runtime/save_shell_state.sh"
)

## Ensure that PASH_TMP_PREFIX is set by pa.sh
assert not os.getenv("PASH_TMP_PREFIX") is None
PASH_TMP_PREFIX = os.getenv("PASH_TMP_PREFIX")

SOCKET_BUF_SIZE = 8192


##
## Global configuration used by all pash components
##
LOGGING_PREFIX = ""
OUTPUT_TIME = False
DEBUG_LEVEL = 0
LOG_FILE = ""


HDFS_PREFIX = "$HDFS_DATANODE_DIR/"


config = {}
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

    ## Also set logging here
    # Format logging
    # ref: https://docs.python.org/3/library/logging.html#formatter-objects
    ## TODO: When we add more logging levels bring back the levelname+time
    if given_pash_args.log_file == "":
        logging.basicConfig(format="%(message)s")
    else:
        logging.basicConfig(
            format="%(message)s",
            filename=f"{os.path.abspath(given_pash_args.log_file)}",
            filemode="w",
        )

    # Set debug level
    if given_pash_args.debug == 1:
        logging.getLogger().setLevel(logging.INFO)
    elif given_pash_args.debug >= 2:
        logging.getLogger().setLevel(logging.DEBUG)


## Increase the recursion limit (it seems that the parser/unparser needs it for bigger graphs)
sys.setrecursionlimit(10000)


def load_config(config_file_path=""):
    global config
    pash_config = {}
    CONFIG_KEY = "distr_planner"

    if config_file_path == "":
        config_file_path = "{}/compiler/config.json".format(PASH_TOP)
    with open(config_file_path) as config_file:
        pash_config = json.load(config_file)

    if not pash_config:
        raise Exception(
            "No valid configuration could be loaded from {}".format(config_file_path)
        )

    if CONFIG_KEY not in pash_config:
        raise Exception(
            "Missing `{}` config in {}".format(CONFIG_KEY, config_file_path)
        )

    config = pash_config


def getWidth():
    cpus = os.cpu_count()
    return math.floor(cpus / 8) if cpus >= 16 else 2


def add_general_config_arguments(parser):
    ## TODO: Delete that at some point, or make it have a different use (e.g., outputting time even without -d 1).
    parser.add_argument(
        "-t",
        "--output_time",  # FIXME: --time
        help="(obsolete, time is always logged now) output the time it took for every step",
        action="store_true",
    )
    parser.add_argument(
        "-d",
        "--debug",
        type=int,
        help="configure debug level; defaults to 0",
        default=0,
    )
    parser.add_argument(
        "--log_file",
        help="configure where to write the log; defaults to stderr.",
        default="",
    )


## These are arguments that are common to pash.py and pash_compiler.py
def add_common_arguments(parser):
    add_general_config_arguments(parser)

    parser.add_argument(
        "-w",
        "--width",
        type=int,
        default=getWidth(),
        help="set data-parallelism factor",
    )
    parser.add_argument(
        "--no_optimize",
        help="not apply transformations over the DFG",
        action="store_true",
    )
    parser.add_argument(
        "--dry_run_compiler",
        help="not execute the compiled script, even if the compiler succeeded",
        action="store_true",
    )
    parser.add_argument(
        "--assert_compiler_success",
        help="assert that the compiler succeeded (used to make tests more robust)",
        action="store_true",
    )
    parser.add_argument(
        "--avoid_pash_runtime_completion",
        help="avoid the pash_runtime execution completion (only relevant when --debug > 0)",
        action="store_true",
    )
    parser.add_argument(
        "--profile_driven",
        help="(experimental) use profiling information when optimizing",
        action="store_true",
    )
    parser.add_argument(
        "-p",
        "--output_optimized",  # FIXME: --print
        help="output the parallel shell script for inspection",
        action="store_true",
    )
    parser.add_argument(
        "--graphviz",
        help="generates graphical representations of the dataflow graphs. The option argument corresponds to the format. PaSh stores them in a timestamped directory in the argument of --graphviz_dir",
        choices=["no", "dot", "svg", "pdf", "png"],
        default="no",
    )
    ## TODO: To discuss: Do we maybe want to have graphviz to always be included
    ##       in the temp directory (under a graphviz subdirectory) instead of in its own?
    ##   kk: I think that ideally we want a log-directory where we can put logs, graphviz,
    ##       and other observability and monitoring info (instead of putting them in the temp).
    parser.add_argument(
        "--graphviz_dir",
        help="the directory in which to store graphical representations",
        default="/tmp",
    )
    parser.add_argument(
        "--no_eager",
        help="(experimental) disable eager nodes before merging nodes",
        action="store_true",
    )
    parser.add_argument(
        "--no_daemon",
        help="(obsolete) does nothing -- Run the compiler everytime we need a compilation instead of using the daemon",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--parallel_pipelines",
        help="(obsolete) Run multiple pipelines in parallel if they are safe to run. Now true by default. See --no_parallel_pipelines.",
        action="store_true",
        default=True,
    )
    parser.add_argument(
        "--no_parallel_pipelines",
        help="Disable parallel running of independent pipelines",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--parallel_pipelines_limit",
        help="Maximum number of parallel independent pipelines",
        type=int,
        default=2,
    )
    parser.add_argument(
        "--r_split_batch_size",
        type=int,
        help="configure the batch size of r_split (default: 1MB)",
        default=1000000,
    )
    parser.add_argument(
        "--r_split",
        help="(obsolete) does nothing -- only here for old interfaces (not used anywhere in the code)",
        action="store_true",
    )
    parser.add_argument(
        "--dgsh_tee",
        help="(obsolete) does nothing -- only here for old interfaces (not used anywhere in the code)",
        action="store_true",
    )
    parser.add_argument(
        "--speculative",
        help="(experimental) use the speculative execution preprocessing and runtime (NOTE: this has nothing to do with --speculation, which is actually misnamed, and should be named concurrent compilation/execution and is now obsolete)",
        action="store_true",
        default=False,
    )
    ## This is misnamed, it should be named concurrent compilation/execution
    parser.add_argument(
        "--speculation",
        help="(obsolete) does nothing -- run the original script during compilation; if compilation succeeds, abort the original and run only the parallel (quick_abort) (Default: no_spec)",
        choices=["no_spec", "quick_abort"],
        default="no_spec",
    )
    parser.add_argument(
        "--termination",
        help="(experimental) determine the termination behavior of the DFG. Defaults to cleanup after the last process dies, but can drain all streams until depletion",
        choices=["clean_up_graph", "drain_stream"],
        default="clean_up_graph",
    )
    parser.add_argument(
        "--daemon_communicates_through_unix_pipes",
        help="(experimental) the daemon communicates through unix pipes instead of sockets",
        action="store_true",
    )
    parser.add_argument(
        "--distributed_exec",
        help="(experimental) execute the script in a distributed environment. Remote machines should be configured and ready",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--config_path",
        help="determines the config file path. By default it is 'PASH_TOP/compiler/config.yaml'.",
        default="",
    )
    parser.add_argument(
        "--version",
        action="version",
        version="%(prog)s {version}".format(version=__version__),
    )
    return


def pass_common_arguments(pash_arguments):
    arguments = []
    if pash_arguments.no_optimize:
        arguments.append("--no_optimize")
    if pash_arguments.dry_run_compiler:
        arguments.append("--dry_run_compiler")
    if pash_arguments.assert_compiler_success:
        arguments.append("--assert_compiler_success")
    if pash_arguments.avoid_pash_runtime_completion:
        arguments.append("--avoid_pash_runtime_completion")
    if pash_arguments.profile_driven:
        arguments.append("--profile_driven")
    if pash_arguments.output_time:
        arguments.append("--output_time")
    if pash_arguments.output_optimized:
        arguments.append("--output_optimized")
    arguments.append("--graphviz")
    arguments.append(pash_arguments.graphviz)
    arguments.append("--graphviz_dir")
    arguments.append(pash_arguments.graphviz_dir)
    if not pash_arguments.log_file == "":
        arguments.append("--log_file")
        arguments.append(pash_arguments.log_file)
    if pash_arguments.no_eager:
        arguments.append("--no_eager")
    if pash_arguments.distributed_exec:
        arguments.append("--distributed_exec")
    if pash_arguments.speculative:
        arguments.append("--speculative")
    if pash_arguments.no_parallel_pipelines:
        arguments.append("--no_parallel_pipelines")
    if pash_arguments.daemon_communicates_through_unix_pipes:
        arguments.append("--daemon_communicates_through_unix_pipes")
    arguments.append("--parallel_pipelines_limit")
    arguments.append(str(pash_arguments.parallel_pipelines_limit))
    arguments.append("--r_split_batch_size")
    arguments.append(str(pash_arguments.r_split_batch_size))
    arguments.append("--debug")
    arguments.append(str(pash_arguments.debug))
    arguments.append("--termination")
    arguments.append(pash_arguments.termination)
    arguments.append("--width")
    arguments.append(str(pash_arguments.width))
    if not pash_arguments.config_path == "":
        arguments.append("--config_path")
        arguments.append(pash_arguments.config_path)
    return arguments


def init_log_file():
    global LOG_FILE
    if not LOG_FILE == "":
        with open(LOG_FILE, "w") as f:
            pass


##
## Set the shell variables
##


def set_vars_file(var_file_path: str, var_dict: dict):
    global config
    config["shell_variables"] = var_dict
    config["shell_variables_file_path"] = var_file_path
