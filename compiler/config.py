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

config = {}

def load_config():
    global config
    dish_config = {}
    CONFIG_KEY = 'distr_planner'

    ## TODO: allow this to be passed as an argument
    config_file_path = '{}/compiler/config.yaml'.format(DISH_TOP)
    with open(config_file_path) as config_file:
        dish_config = yaml.load(config_file, Loader=yaml.FullLoader)

    if not dish_config:
        raise Exception('No valid configuration could be loaded from {}'.format(config_file_path))

    if CONFIG_KEY not in dish_config:
        raise Exception('Missing `{}` config in {}'.format(CONFIG_KEY, config_file_path))

    config = dish_config[CONFIG_KEY]
