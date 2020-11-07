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
    compilation_start_time = datetime.now()
    ## Parse arguments
    args = parse_args()
    config.pash_args = args

    if not config.config:
        config.load_config(args.config_path)

    ## Load annotations
    ## TODO: This should not happen here anymore.
    config.load_annotation_files(config.config['distr_planner']['annotations_dir'])

    ## 1. Execute the POSIX shell parser that returns the AST in JSON
    input_script_path = args.input
    json_ast_string = parse_shell(input_script_path)

    ## 2. Parse JSON to AST objects
    ast_objects = parse_json_ast_string(json_ast_string)

    ## 3. Preprocess ASTs by replacing possible candidates for compilation
    ##    with calls to the PaSh runtime.
    compiled_asts = preprocess(ast_objects, config.config)

    ## 4. Translate the new AST back to shell syntax
    input_script_wo_extension, input_script_extension = os.path.splitext(input_script_path)
    ir_filename = input_script_wo_extension + ".ir"
    save_asts_json(compiled_asts, ir_filename)
    from_ir_to_shell_file(ir_filename, args.output)

    compilation_end_time = datetime.now()
    print_time_delta("Preprocessing", compilation_start_time, compilation_end_time, args)

    ## TODO: Change all occurences of compile to preprocess since PaSh just 
    #        preprocesses and replaces possibly parallelizable regions with calls to PaSh.

    ## 5. Execute the preprocessed version of the input script
    if(not args.compile_only):
        execute_script(args.output, args.output_optimized, args.compile_optimize_only)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="the script to be compiled and executed")
    parser.add_argument("output", help="the path of the compiled shell script")
    parser.add_argument("--compile_only", help="only compile the input script and not execute it",
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

    ## TODO: Delete when done
    # final_asts = []
    # ## This is for the files in the IR
    # fileIdGen = FileIdGen()

    # ## Compile the asts
    # compiled_asts = compile_asts(ast_objects, fileIdGen, config)

    # for i, compiled_ast in enumerate(compiled_asts):
    #     # print("Replacing AST {}".format(i))
    #     # print(compiled_ast)

    #     ## Replace the IRs in the ASTs with calls to the distribution
    #     ## planner. Save the IRs in temporary files.
    #     final_ast = replace_irs(compiled_ast, irFileGen, config)

    #     # print("Final AST:")
    #     # print(final_ast)
    #     final_asts.append(final_ast)
    # return final_asts

def execute_script(compiled_script_filename, output_optimized, compile_optimize_only):
    exec_obj = subprocess.run(["/bin/bash", compiled_script_filename])
    exec_obj.check_returncode()

if __name__ == "__main__":
    main()
