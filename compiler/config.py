import os
import subprocess
import yaml

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

    config = dish_config[CONFIG_KEY]


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
