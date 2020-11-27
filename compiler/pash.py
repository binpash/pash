import os
import subprocess
import argparse
from datetime import datetime

from ast_to_ir import *
from pash_runtime import *
from ir import *
from json_ast import *
from parse import parse_shell, from_ir_to_shell_file
from util import *
import config

def main():
    preprocessing_start_time = datetime.now()
    ## Parse arguments
    args = parse_args()
    config.pash_args = args

    if not config.config:
        config.load_config(args.config_path)

    ## Load annotations
    ## TODO: The annotations are not used in the preprocessing step anymore
    ##       so we can avoid loading the annotation files here.
    # config.load_annotation_files(config.config['distr_planner']['annotations_dir'])

    ## 1. Execute the POSIX shell parser that returns the AST in JSON
    input_script_path = args.input
    json_ast_string = parse_shell(input_script_path)

    ## 2. Parse JSON to AST objects
    ast_objects = parse_json_ast_string(json_ast_string)

    ## 3. Preprocess ASTs by replacing possible candidates for compilation
    ##    with calls to the PaSh runtime.
    preprocessed_asts = preprocess(ast_objects, config.config)

    ## 4. Translate the new AST back to shell syntax
    input_script_wo_extension, _input_script_extension = os.path.splitext(input_script_path)
    ir_filename = input_script_wo_extension + ".ir"
    save_asts_json(preprocessed_asts, ir_filename)
    from_ir_to_shell_file(ir_filename, args.output)

    preprocessing_end_time = datetime.now()
    print_time_delta("Preprocessing", preprocessing_start_time, preprocessing_end_time, args)

    ## 5. Execute the preprocessed version of the input script
    if(not args.preprocess_only):
        execute_script(args.output)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="the script to be compiled and executed")
    parser.add_argument("output", help="the path of the compiled shell script")
    parser.add_argument("--preprocess_only", help="only preprocess the input script and not execute it",
                        action="store_true")
    config.add_common_arguments(parser)
    args = parser.parse_args()
    return args

def preprocess(ast_objects, config):
    ## This is ids for the remporary files that we will save the IRs in
    irFileGen = FileIdGen()

    ## Preprocess ASTs by replacing AST regions with calls to PaSh's runtime. 
    ## Then the runtime will do the compilation and optimization with additional 
    ## information.
    preprocessed_asts = replace_ast_regions(ast_objects, irFileGen, config)

    return preprocessed_asts

def execute_script(compiled_script_filename):
    exec_obj = subprocess.run(["/bin/bash", compiled_script_filename])
    exec_obj.check_returncode()

if __name__ == "__main__":
    main()
