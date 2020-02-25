import os
import subprocess
import argparse

from ast_to_ir import *
from distr_plan import *
from ir import *
from json_ast import *
import config

def main():
    ## Parse arguments
    args = parse_args()

    # global config
    if not config.config:
        config.load_config()

    ## 1. Execute the POSIX shell parser that returns the AST in JSON
    input_script_path = args.input
    json_ast_string = parse_shell(input_script_path)

    ## 2. Parse JSON to AST objects
    ast_objects = parse_json_ast_string(json_ast_string)

    ## 3. Compile ASTs to our intermediate representation
    compiled_asts = compile_asts(ast_objects, config.config)

    ## 4. Translate the new AST back to shell syntax
    input_script_wo_extension, input_script_extension = os.path.splitext(input_script_path)
    ir_filename = input_script_wo_extension + ".ir"
    save_asts_json(compiled_asts, ir_filename)
    from_ir_to_shell(ir_filename, args.output)

    ## 5. Execute the compiled version of the input script
    if(not args.compile_only):
        execute_script(args.output, args.output_optimized, args.compile_optimize_only)

def parse_shell(input_script_path):
    parser_output = subprocess.run([config.PARSER_BINARY, input_script_path], capture_output=True, text=True)
    parser_output.check_returncode()
    return parser_output.stdout

def from_ir_to_shell(ir_filename, new_shell_filename):
    printer_output = subprocess.run([config.PRINTER_BINARY, ir_filename], capture_output=True, text=True)
    printer_output.check_returncode()
    compiled_script = printer_output.stdout
    with open(new_shell_filename, 'w') as new_shell_file:
        new_shell_file.write(compiled_script)

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="the script to be compiled and executed")
    parser.add_argument("output", help="the path of the compiled shell script")
    parser.add_argument("--compile_only", help="only compile the input script and not execute it",
                        action="store_true")
    parser.add_argument("--compile_optimize_only", help="only compile and optimize the input script and not execute it",
                        action="store_true")
    parser.add_argument("--output_optimized", help="output the optimized shell script that was produced by the planner for inspection",
                        action="store_true")
    args = parser.parse_args()
    return args

def compile_asts(ast_objects, config):
    ## This is for the files in the IR
    fileIdGen = FileIdGen()

    ## This is ids for the remporary files that we will save the IRs in
    irFileGen = FileIdGen()

    final_asts = []
    for i, ast_object in enumerate(ast_objects):
        # print("Compiling AST {}".format(i))
        # print(ast_object)

        ## Compile subtrees of the AST to out intermediate representation
        compiled_ast = compile_ast(ast_object, fileIdGen, config)

        # print("Compiled AST:")
        # print(compiled_ast)

        ## Replace the IRs in the ASTs with calls to the distribution
        ## planner. Save the IRs in temporary files.
        final_ast = replace_irs(compiled_ast, irFileGen, config)

        # print("Final AST:")
        # print(final_ast)
        final_asts.append(final_ast)
    return final_asts

def execute_script(compiled_script_filename, output_optimized, compile_optimize_only):
    exec_obj = subprocess.run(["/bin/bash", compiled_script_filename])
    exec_obj.check_returncode()

if __name__ == "__main__":
    main()
