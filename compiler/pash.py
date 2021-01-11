import os
import subprocess
import argparse
from datetime import datetime

from annotations import *
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
    
    ## Initialize the log file
    config.init_log_file()

    if not config.config:
        config.load_config(args.config_path)

    ## Load annotations
    ## TODO: The annotations are not used in the preprocessing step anymore
    ##       so we can avoid loading the annotation files here.
    # config.annotations = load_annotation_files(config.config['distr_planner']['annotations_dir'])

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
    input_script_basename = os.path.basename(input_script_wo_extension)
    ir_filename = os.path.join("/tmp", get_pash_prefixed_random_string() + "_" + input_script_basename + ".ir")
    save_asts_json(preprocessed_asts, ir_filename)

    preprocessed_output_filename = os.path.join("/tmp", get_pash_prefixed_random_string())
    log("Preprocessed script stored in:", preprocessed_output_filename)
    if(args.output_preprocessed):
        log("Preprocessed script:")
        log(from_ir_to_shell(ir_filename))
    from_ir_to_shell_file(ir_filename, preprocessed_output_filename)

    preprocessing_end_time = datetime.now()
    print_time_delta("Preprocessing", preprocessing_start_time, preprocessing_end_time, args)

    ## 5. Execute the preprocessed version of the input script
    if(not args.preprocess_only):
        execute_script(preprocessed_output_filename)


def parse_args():
    prog_name = sys.argv[0]
    if 'PASH_FROM_SH' in os.environ:
        prog_name = os.environ['PASH_FROM_SH']
    parser = argparse.ArgumentParser(prog_name)
    parser.add_argument("input", help="the script to be compiled and executed")
    parser.add_argument("--preprocess_only", 
                        help="only preprocess the input script and not execute it",
                        action="store_true")
    parser.add_argument("--output_preprocessed", 
                        help=" output the preprocessed script",
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
    ## Return the exit code of the executed script
    exit(exec_obj.returncode)

if __name__ == "__main__":
    main()
