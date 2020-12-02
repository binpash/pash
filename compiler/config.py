import json
import os
import subprocess
import yaml

from ir_utils import *

## Global
GIT_TOP_CMD = [ 'git', 'rev-parse', '--show-toplevel', '--show-superproject-working-tree']
if 'PASH_TOP' in os.environ:
    PASH_TOP = os.environ['PASH_TOP']
else:
    PASH_TOP = subprocess.run(GIT_TOP_CMD, capture_output=True,
                              text=True).stdout.rstrip()

PARSER_BINARY = os.path.join(PASH_TOP, "parser/parse_to_json.native")
PRINTER_BINARY = os.path.join(PASH_TOP, "parser/json_to_shell.native")

PYTHON_VERSION = "python3.8"
PLANNER_EXECUTABLE = os.path.join(PASH_TOP, "compiler/pash_runtime.py")
RUNTIME_EXECUTABLE = os.path.join(PASH_TOP, "compiler/pash_runtime.sh")

config = {}
annotations = []
pash_args = None

def load_config(config_file_path=""):
    global config
    pash_config = {}
    CONFIG_KEY = 'distr_planner'

    if(config_file_path == ""):
      config_file_path = '{}/compiler/config.yaml'.format(PASH_TOP)
    with open(config_file_path) as config_file:
        pash_config = yaml.load(config_file, Loader=yaml.FullLoader)

    if not pash_config:
        raise Exception('No valid configuration could be loaded from {}'.format(config_file_path))

    if CONFIG_KEY not in pash_config:
        raise Exception('Missing `{}` config in {}'.format(CONFIG_KEY, config_file_path))

    config = pash_config

## These are arguments that are common to pash.py and pash_runtime.py
def add_common_arguments(parser):
    parser.add_argument("--compile_only", help="only preprocess and compile the input script and not execute it",
                        action="store_true")
    parser.add_argument("--compile_optimize_only",
                        help="only preprocess, compile, and optimize the input script and not execute it",
                        action="store_true")
    parser.add_argument("--output_time", help="output the the time it took for every step in stderr",
                        action="store_true")
    parser.add_argument("--output_optimized",
                        help="output the optimized shell script that"
                        "was produced by the planner for inspection",
                        action="store_true")
    parser.add_argument("--no_eager",
                        help="disable eager nodes before merging nodes",
                        action="store_true")
    parser.add_argument("--termination",
                        help="determines the termination behavior of the DFG. By default it cleans up the graph after the final node dies but it can also drain all streams until depletion.",
                        choices=['clean_up_graph', 'drain_stream'],
                        default="clean_up_graph")
    parser.add_argument("--split_fan_out",
                        type=int,
                        default=1,
                        help="determines the fan out of inserted splits in the DFG")
    parser.add_argument("--config_path",
                        help="determines the config file path. By default it is 'PASH_TOP/compiler/config.yaml'.",
                        default="")
    return

def pass_common_arguments(pash_arguments):
    arguments = []
    if (pash_arguments.compile_only):
        arguments.append(string_to_argument("--compile_only"))
    if (pash_arguments.compile_optimize_only):
        arguments.append(string_to_argument("--compile_optimize_only"))
    if (pash_arguments.output_time):
        arguments.append(string_to_argument("--output_time"))
    if (pash_arguments.output_optimized):
        arguments.append(string_to_argument("--output_optimized"))
    if (pash_arguments.no_eager):
        arguments.append(string_to_argument("--no_eager"))
    arguments.append(string_to_argument("--termination"))
    arguments.append(string_to_argument(pash_arguments.termination))
    arguments.append(string_to_argument("--split_fan_out"))
    arguments.append(string_to_argument(str(pash_arguments.split_fan_out)))
    if(not pash_arguments.config_path == ""):
        arguments.append(string_to_argument("--config_path"))
        arguments.append(string_to_argument(pash_arguments.config_path))
    return arguments

##
## Load annotation files
##

def load_annotation_file(abs_annotation_filename):
    with open(abs_annotation_filename) as annotation_file:
        try:
            annotation = json.load(annotation_file)
            return [annotation]
        except json.JSONDecodeError as err:
            print("WARNING: Could not parse annotation for file:", abs_annotation_filename)
            print("|-- {}".format(err))
            return []

def load_annotation_files(annotation_dir):
    global annotations
    if(not os.path.isabs(annotation_dir)):
        annotation_dir = os.path.join(PASH_TOP, annotation_dir)

    for (dirpath, dirnames, filenames) in os.walk(annotation_dir):
        json_filenames = [os.path.join(dirpath, filename) for filename in filenames
                          if filename.endswith(".json")]
        curr_annotations = [ann for filename in json_filenames for ann in load_annotation_file(filename) ]
        annotations += curr_annotations


# TODO load command class file path from config
command_classes_file_path = '{}/compiler/command-classes.yaml'.format(PASH_TOP)
command_classes = {}
with open(command_classes_file_path) as command_classes_file:
    command_classes = yaml.load(command_classes_file, Loader=yaml.FullLoader)

if not command_classes:
    raise Exception('Failed to load description of command classes from {}'.format(command_classes_file_path))

stateless_commands = command_classes['stateless'] if 'stateless' in command_classes else {}
pure_commands = command_classes['pure'] if 'pure' in command_classes else {}
parallelizable_pure_commands = command_classes['parallelizable_pure'] if 'parallelizable_pure' in command_classes else {}

## TODO: Move these to a configuration file
bigram_g_map_num_outputs = 3
