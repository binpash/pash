import argparse
import traceback

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

def parse_command(input):
    ## TODO: Improve the way parsing happens plz :')
    if(input.startswith("Compile:")):
        return compile(input)
    else:
        return error_response(f'Unsupported command: {input}')

## TODO: Improve the way parsing happens plz :') At the moment this will fail with : in file etc
def parse_compile_command(input):
    try:
        components = input.rstrip().split("|")
        compiled_script_file = components[0].split(":")[1]
        var_file = components[1].split(":")[1]
        input_ir_file = components[2].split(":")[1]
        return compiled_script_file, var_file, input_ir_file
    except:
        raise Exception(f'Parsing failure for line: {input}')

def compile(input):
    compiled_script_file, var_file, input_ir_file = parse_compile_command(input)
    
    ## Read any shell variables files if present
    config.read_vars_file(var_file)

    ## Call the main procedure
    pash_runtime.compile_optimize_output_script(input_ir_file, compiled_script_file, config.pash_args)

    return success_response(f'{compiled_script_file} {var_file} {input_ir_file}')

def main():
    args = init()

    while True:
        ## Process a single request
        with open(args.input) as fin, open(args.output, "w") as fout:
            input = fin.read()
            try:
                ret = parse_command(input)
            except Exception as e:
                log(traceback.format_exc())
                ret = error_response(str(e))
            
            fout.write(ret)
            fout.flush()

if __name__ == "__main__":
    main()
