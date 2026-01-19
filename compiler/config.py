import json
import logging
import os
import subprocess

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

BASH_VERSION = tuple(int(i) for i in os.getenv("PASH_BASH_VERSION").split(" "))


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
