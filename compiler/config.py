import json
import os
import subprocess
import yaml

from ir_utils import *

## Global
GIT_TOP_CMD = [ 'git', 'rev-parse', '--show-toplevel', '--show-superproject-working-tree']
if 'DISH_TOP' in os.environ:
    DISH_TOP = os.environ['DISH_TOP']
else:
    DISH_TOP = subprocess.run(GIT_TOP_CMD, capture_output=True,
            text=True).stdout.rstrip()

PARSER_BINARY = os.path.join(DISH_TOP, "parser/parse_to_json.native")
PRINTER_BINARY = os.path.join(DISH_TOP, "parser/json_to_shell.native")

PYTHON_VERSION = "python3.8"
PLANNER_EXECUTABLE = os.path.join(DISH_TOP, "compiler/distr_plan.py")

config = {}
annotations = []
dish_args = None

def load_config(config_file_path=False):
    global config
    dish_config = {}
    CONFIG_KEY = 'distr_planner'

    if not config_file_path:
      config_file_path = '{}/compiler/config.yaml'.format(DISH_TOP)
    with open(config_file_path) as config_file:
        dish_config = yaml.load(config_file, Loader=yaml.FullLoader)

    if not dish_config:
        raise Exception('No valid configuration could be loaded from {}'.format(config_file_path))

    if CONFIG_KEY not in dish_config:
        raise Exception('Missing `{}` config in {}'.format(CONFIG_KEY, config_file_path))

    config = dish_config

## These are arguments that are common to dish and distr_plan.py
def add_common_arguments(parser):
    parser.add_argument("--compile_optimize_only",
                        help="only compile and optimize the input script and not execute it",
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
    parser.add_argument("--clean_up_graph",
                        help="clean up the parallel dataflow graphs when the final node dies",
                        action="store_true")
    parser.add_argument("--drain_streams",
                        help="drain up all streams instead of letting them hang (this is an alternative to --clean_up_graph",
                        action="store_true")
    parser.add_argument("--auto_split",
                        help="uses a no-task-parallelism split that automatically calculates batch size",
                        action="store_true")
    return

def pass_common_arguments(dish_arguments):
    arguments = []
    if (dish_arguments.compile_optimize_only):
        arguments.append(string_to_argument("--compile_optimize_only"))
    if (dish_arguments.output_time):
        arguments.append(string_to_argument("--output_time"))
    if (dish_arguments.output_optimized):
        arguments.append(string_to_argument("--output_optimized"))
    if (dish_arguments.no_eager):
        arguments.append(string_to_argument("--no_eager"))
    if (dish_arguments.clean_up_graph):
        arguments.append(string_to_argument("--clean_up_graph"))
    if (dish_arguments.drain_streams):
        arguments.append(string_to_argument("--drain_streams"))
    if (dish_arguments.auto_split):
        arguments.append(string_to_argument("--auto_split"))
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
        annotation_dir = os.path.join(DISH_TOP, annotation_dir)

    for (dirpath, dirnames, filenames) in os.walk(annotation_dir):
        json_filenames = [os.path.join(dirpath, filename) for filename in filenames
                          if filename.endswith(".json")]
        curr_annotations = [ann for filename in json_filenames for ann in load_annotation_file(filename) ]
        annotations += curr_annotations


# TODO load command class file path from config
command_classes_file_path = '{}/compiler/command-classes.yaml'.format(DISH_TOP)
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
