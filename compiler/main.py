import os
import subprocess
import argparse

from ast_to_ir import *
from distr_plan import *
from ir import *
from json_ast import *

PARSER_VAR = "DISH_PARSER"

def main():
    ## Translation process:
    ##
    ## 1. Parse json to an AST object
    ##
    ## 2. TODO: Ensure that the AST is in our supported subset (TODO)
    ##
    ## 3. Compile subtrees of the AST to out intermediate representation
    ##
    ## 4. Replace the IRs in the ASTs with calls to the distribution
    ##    planner. Save the IRs in temporary files.
    ##
    ## 5. Translate the new AST back to shell syntax
    ##
    ## 6. Execute the new AST using a standard shell

    args = parse_args()

    input_script_path = args.input
    parser_binary = os.environ[PARSER_VAR]

    ## Execute the POSIX shell parser that returns the AST in JSON
    parser_output = subprocess.run([parser_binary, input_script_path], capture_output=True, text=True)
    parser_output.check_returncode()

    json_ast_string = parser_output.stdout

    # ast_objects = parse_json_ast(json_filename)
    ast_objects = parse_json_ast_string(json_ast_string)
    check_if_asts_supported(ast_objects)

    ## This is for the files in the IR
    fileIdGen = FileIdGen()

    ## This is ids for the remporary files that we will save the IRs in
    irFileGen = FileIdGen()

    final_asts = []
    for i, ast_object in enumerate(ast_objects):
        # print("Compiling AST {}".format(i))
        # print(ast_object)
        compiled_ast = compile_ast(ast_object, fileIdGen)
        # print("Compiled AST:")
        # print(compiled_ast)

        final_ast = replace_irs(compiled_ast, irFileGen)
        # print("Final AST:")
        # print(final_ast)
        final_asts.append(final_ast)

    ## TODO: The following lines are currently useless, since we just
    ## execute the first dataflow in each script manually
    ir_filename = input_script_path + ".ir"
    save_asts_json(final_asts, ir_filename)
    new_shell_filename = input_script_path + "_compiled.sh"
    from_ir_to_shell(ir_filename, new_shell_filename)

    ## Execute the compiled version of the input script
    if(not args.compile_only):
        execute_script(new_shell_filename, args.output)

## Note: It seems that in order for distribution and planning to
## happen correctly, the planner has to be invoked as late as possible
## during the shell execution.

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
    args = parser.parse_args()
    return args

## TODO: Extend this to properly execute the compiled script
def execute_script(compiled_script_filename, output_script_path):
    ## For now, just execute the first dataflow graph in the script
    ir_filename = "/tmp/dish_temp_ir_file1"
    output_dir = "/tmp/"
    fan_out = 4
    batch_size = 100
    optimize_script(ir_filename, output_script_path, output_dir, fan_out, batch_size)

##
## Ideally, we would like to execute part of the shell script
## normally, and every time we reach the IR, we would like to pass
## control to our distribution planner. The distribution planner might
## need to ask some queries to the shell and DFS, like values of
## variables, file locations, etc..., and then it should be able to
## plan how to distribute the query. So essentially, we would like to
## replace every IR node in the AST with a command (or function call)
## that gathers as much information as possible from the shell state,
## and then calls the planner (could be a simple call like `python3
## distr_plan.py ...`). The planner then plans how to distribute an IR
## in different nodes (after having received in its arguments
## locations of files etc).
##

if __name__ == "__main__":
    main()
