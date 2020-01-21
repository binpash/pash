import os
import subprocess
import argparse

from ast_to_ir import *
from distr_plan import *
from ir import *
from json_ast import *

GIT_TOP_CMD = [ 'git', 'rev-parse', '--show-toplevel', '--show-superproject-working-tree']
DISH_TOP = os.environ['DISH_TOP'] or subprocess.run(GIT_TOP_CMD, capture_output=True, text=True).stdout
PARSER_VAR = "DISH_PARSER" # TODO (nv): could pick from DISH_TOP

def main():
    ## Parse arguments
    args = parse_args()

    ## 1. Execute the POSIX shell parser that returns the AST in JSON
    input_script_path = args.input
    parser_binary = os.environ[PARSER_VAR]
    parser_output = subprocess.run([parser_binary, input_script_path], capture_output=True, text=True)
    parser_output.check_returncode()
    json_ast_string = parser_output.stdout

    ## 2. Parse JSON to AST objects
    ast_objects = parse_json_ast_string(json_ast_string)

    ## 3. Compile ASTs to our intermediate representation
    compiled_asts = compile_asts(ast_objects)

    ## TODO: The following lines are currently useless, since we just
    ## execute the first dataflow in each script manually
    ## TODO: Don't hardcode the .ir file name
    ## 4. TODO: Translate the new AST back to shell syntax
    ir_filename = input_script_path + ".ir"
    save_asts_json(compiled_asts, ir_filename)
    new_shell_filename = input_script_path + "_compiled.sh"
    from_ir_to_shell(ir_filename, new_shell_filename)

    ## 5. Execute the compiled version of the input script
    if(not args.compile_only):
        execute_script(new_shell_filename, args.output, args.output_optimized)

def from_ir_to_shell(ir_filename, new_shell_filename):
    ## TODO: Execute the ocaml json_to_shell.native and save its
    ## output in a file
    return

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="the script to be compiled and executed")
    parser.add_argument("output", help="the path of the compiled shell script")
    parser.add_argument("--compile_only", help="only compile the input script and not execute it",
                        action="store_true")
    parser.add_argument("--output_optimized", help="output the optimized shell script that was produced by the planner for inspection",
                        action="store_true")
    args = parser.parse_args()
    return args

def compile_asts(ast_objects):
    ## This is for the files in the IR
    fileIdGen = FileIdGen()

    ## This is ids for the remporary files that we will save the IRs in
    irFileGen = FileIdGen()

    final_asts = []
    for i, ast_object in enumerate(ast_objects):
        # print("Compiling AST {}".format(i))
        # print(ast_object)

        ## Compile subtrees of the AST to out intermediate representation
        compiled_ast = compile_ast(ast_object, fileIdGen)

        # print("Compiled AST:")
        # print(compiled_ast)

        ## Replace the IRs in the ASTs with calls to the distribution
        ## planner. Save the IRs in temporary files.
        final_ast = replace_irs(compiled_ast, irFileGen)

        # print("Final AST:")
        # print(final_ast)
        final_asts.append(final_ast)
    return final_asts

## TODO: Extend this to properly execute the compiled script
def execute_script(compiled_script_filename, output_script_path, output_optimized):
    ## For now, just execute the first dataflow graph in the script
    optimize_script(output_script_path)

if __name__ == "__main__":
    main()
