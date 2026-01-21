import json
import logging
import os
import subprocess
import sys


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
SAVE_ARGS_EXECUTABLE = os.path.join(PASH_TOP, "runtime/save_args.sh")

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
    if given_pash_args.debug == 0:
        logging.getLogger().setLevel(logging.ERROR)
    elif given_pash_args.debug == 1:
        logging.getLogger().setLevel(logging.WARNING)
    elif given_pash_args.debug == 2:
        logging.getLogger().setLevel(logging.INFO)
    elif given_pash_args.debug >= 3:
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
