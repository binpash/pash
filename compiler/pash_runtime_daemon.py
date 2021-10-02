import argparse

from annotations import *
import config
import pash_runtime
from util import *

##
## A Daemon responding to requests for compilation
##
## Note: Not an actual daemon with the strict Unix sense
##

## TODO: Rename the pash_runtime to pash_compiler and this to pash_daemon

## TODO: Should we maybe use sockets instead of fifos?

## TODO: Fix the daemon logging.

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="the input fifo from which the daemon will read its input")
    parser.add_argument("output", help="the output fifo to which the daemon will write its output")
    config.add_common_arguments(parser)
    args, unknown_args = parser.parse_known_args()

    ## Print all the arguments before they are modified below
    # print("Arguments:")
    # for arg_name, arg_val in vars(args).items():
    #     print(arg_name, arg_val)
    # print("-" * 40)
    # print("Unknown rest arguments:")
    # print(unknown_args)
    # print("-" * 40)

    return args

## Initialize the daemon
def init():
    args = parse_args()
    config.pash_args = args

    ## Load the configuration
    if not config.config:
        config.load_config(args.config_path)
    
    ## Load annotations
    config.annotations = load_annotation_files(config.config['distr_planner']['annotations_dir'])

    pash_runtime.runtime_config = config.config['distr_planner']

    return args

def success_response(string):
    return f'OK: {string}\n'

def error_response(string):
    return f'ERROR: {string}\n'

def parse_command_line(line):
    ## TODO: Improve the way parsing happens plz :')
    if(line.startswith("Compile:")):
        return compile(line)
    else:
        return error_response(f'Unsupported command: {line}')

## TODO: Improve the way parsing happens plz :') At the moment this will fail with : in file etc
def parse_compile_line(line):
    try:
        components = line.rstrip().split("|")
        compiled_script_file = components[0].split(":")[1]
        var_file = components[1].split(":")[1]
        return compiled_script_file, var_file
    except:
        raise Exception(f'Parsing failure for line: {line}')

def compile(line):
    compiled_script_file, var_file = parse_compile_line(line)
    ## TODO: Load the variable file

    ## TODO: Compile and optimize the df region script
    return success_response(f'{compiled_script_file} {var_file}')

def main():
    args = init()

    while True:
        ## Process a single request
        with open(args.input) as fin, open(args.output, "w") as fout:
            for line in fin.readlines():
                try:
                    ret = parse_command_line(line)
                except Exception as e:
                    ret = error_response(str(e))
                    
                fout.write(ret)

if __name__ == "__main__":
    main()