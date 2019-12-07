import os
import subprocess

from ast_to_ir import *
from distr_plan import *
from ir import *
from json_ast import *

PARSER_VAR = "DISH_PARSER"

def from_ir_to_shell(ir_filename, new_shell_filename):
    ## TODO: Execute the ocaml json_to_shell.native and save its
    ## output in a file
    return


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

    input_script_path = sys.argv[1]
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

    ##
    ## TODO: Make sure that the whole end-to-end prototype runs (calling
    ## ocaml to parse, and make a trivial distribution planner that runs
    ## runs the original ast.



## Note: It seems that in order for distribution and planning to
## happen correctly, the planner has to be invoked as late as possible
## during the shell execution.

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
