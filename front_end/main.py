from ast_to_ir import *
from distr_plan import *
from ir import *

# The following files just contain POSIX (no process substitution)
# pipes and commands
#
# json_filename = "../scripts/json/compile.sh.json"
# json_filename = "../scripts/json/grep.sh.json"
json_filename = "../scripts/json/minimal.sh.json"
# json_filename = "../scripts/json/shortest-scripts.sh.json"
# json_filename = "../scripts/json/spell.sh.json"
# json_filename = "../scripts/json/topn.sh.json"
# json_filename = "../scripts/json/wc.sh.json"
# json_filename = "../scripts/json/wf.sh.json"
# json_filename = "../scripts/json/ngrams.sh.json"
# json_filename = "../scripts/json/diff.sh.json" 

# The following contain And operatos together with pipes and commands
#
# json_filename = "../scripts/json/old_ngrams.sh.json" 


# The following is interesting, since it contains command substitution
# (which is parsed as backticks in the Greenberg parser). The
# backticks seem to mean that whatever is in the backticks will
# execute first, and its output will become a string in place of the
# backticks
# json_filename = "../scripts/json/page-count.sh.json"

# Unidentified
#
# json_filename = "../scripts/json/maximal.sh.json"



## Translation process:
## 1. Parse json to an AST object
## 2. Ensure that the AST is in our supported subset (TODO)
## 3. Compile the AST to out intermediate representation

ast_objects = parse_json_ast(json_filename)
check_if_asts_supported(ast_objects)
fileIdGen = FileIdGen()

for i, ast_object in enumerate(ast_objects):
    print("Compiling AST {}".format(i))
    print(ast_object)
    compiled_ast = compile_ast(ast_object, fileIdGen)
    print("Compiled AST:")
    print(compiled_ast)
    

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
## TODO: (If we are following this plan)
##
## - Change compile_ast to save each internal IR to a temporary file
##   (in some format) and then simply do a command to `python3
##   distr_plan.py IR_file`. (Optimization: Think about what information
##   from the shell state would help in distribution)
##
## - Implement a distr_plan.py file that reads an IR and a network
##   configuration from two files, and then plans how to distribute
##   the IR computation.
##
## - Extend the main.py to rewrite the AST in json, and call the
##   libdash pretty printer, to actually produce a new shell
##   script. Then, run this script.
