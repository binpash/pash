import sys
import os
import subprocess
import argparse
from datetime import datetime
from annotations import *
from ast_to_ir import *
from ir import *
from parse import parse_shell_to_asts, parse_shell_to_asts_interactive, from_ast_objects_to_shell
from pash_graphviz import maybe_init_graphviz_dir
from util import *
import config
import shutil

def main():
    ## Set the logging prefix
    config.LOGGING_PREFIX = "PaSh Preprocessor: "

    ## Parse arguments
    args, shell_name = parse_args()
    ## If it is interactive we need a different execution mode
    ##
    ## The user can also ask for an interactive mode irregardless of whether pash was invoked in interactive mode. 
    if(len(args.input) == 0 or args.interactive):
        interactive(args, shell_name)
    else:
        ## 1. Execute the POSIX shell parser that returns the AST in JSON
        input_script_path = args.input[0]
        input_script_arguments = args.input[1:]
        preprocessing_parsing_start_time = datetime.now()
        ast_objects = parse_shell_to_asts(input_script_path)
        preprocessing_parsing_end_time = datetime.now()
        print_time_delta("Preprocessing -- Parsing", preprocessing_parsing_start_time, preprocessing_parsing_end_time, args)

        ## Preprocess and execute the parsed ASTs
        return_code = preprocess_and_execute_asts(ast_objects, args, input_script_arguments, shell_name)
        
        log("-" * 40) #log end marker
        ## Return the exit code of the executed script
        sys.exit(return_code)

def preprocess_ast(ast_objects, args):
    ## 2. Preprocess ASTs by replacing possible candidates for compilation
    ##    with calls to the PaSh runtime.
    preprocessing_pash_start_time = datetime.now()
    preprocessed_asts = preprocess(ast_objects, config.config)
    preprocessing_pash_end_time = datetime.now()
    print_time_delta("Preprocessing -- PaSh", preprocessing_pash_start_time, preprocessing_pash_end_time, args)

    ## 3. Translate the new AST back to shell syntax
    preprocessing_unparsing_start_time = datetime.now()
    preprocessed_shell_script = from_ast_objects_to_shell(preprocessed_asts)
    if(args.output_preprocessed):
        log("Preprocessed script:")
        log(preprocessed_shell_script)
    
    preprocessing_unparsing_end_time = datetime.now()
    print_time_delta("Preprocessing -- Unparsing", preprocessing_unparsing_start_time, preprocessing_unparsing_end_time, args)
    return preprocessed_shell_script


def preprocess_and_execute_asts(ast_objects, args, input_script_arguments, shell_name):
    preprocessed_shell_script = preprocess_ast(ast_objects, args)
    
    ## Write the new shell script to a file to execute
    _, fname = ptempfile()
    log("Preprocessed script stored in:", fname)
    with open(fname, 'w') as new_shell_file:
        new_shell_file.write(preprocessed_shell_script)


    ## 4. Execute the preprocessed version of the input script
    if(not args.preprocess_only):
        return_code = execute_script(fname, args.command, input_script_arguments, shell_name)
    else:
        return_code = 0

    return return_code

## TODO: Create an interactive pash
def interactive(args, shell_name):
    ## This means that we are interactive
    assert(len(args.input) == 0 or args.interactive)

    if len(args.input) == 0:
        ## Read from stdin
        input_script_path = "-"
        script_args = []
    else:
        input_script_path = args.input[0]
        script_args = args.input[1:]

    ## Spawn a bash shell in interactive mode
    new_env = shell_env(shell_name)
    subprocess_args = bash_prefix_args()
    ## Add this option to read code from stdin
    subprocess_args.append("-s")
    ## Add the script arguments here
    subprocess_args += script_args
    with subprocess.Popen(subprocess_args,
                          env=new_env,
                          stdin=subprocess.PIPE,
                          universal_newlines=True,
                          close_fds=False) as shell_proc:
        ## TODO: Do we need to pipe stdout/stderror

        ## First send an exec so that we change the name of the shell
        ##
        ## TODO: Can this be done in a less ad-hoc way?
        command = bash_exec_string(shell_name)
        shell_proc.stdin.write(command)

        ## For each parsed AST:
        ##   1. Preprocess it
        ##   2. Translate it to shell syntax
        ##   3. Send it to the interactive bash 
        ast_objects = parse_shell_to_asts_interactive(input_script_path)
        for ast_object in ast_objects:
            ## Preprocess each ast object and produce a preprocessed shell script fragment
            preprocessed_shell_script = preprocess_ast([ast_object], args)
            log("Sending script to shell process...")
            ## Send the preprocessed script fragment to the shell process
            shell_proc.stdin.write(preprocessed_shell_script)
            shell_proc.stdin.flush()
    
        ## Close the input and wait for the internal process to finish
        shell_proc.stdin.close()
        shell_proc.wait()
        
        log("-" * 40) #log end marker
        ## Return the exit code of the executed script
        sys.exit(shell_proc.returncode)


def parse_args():
    prog_name = sys.argv[0]
    if 'PASH_FROM_SH' in os.environ:
        prog_name = os.environ['PASH_FROM_SH']
    ## We need to set `+` as a prefix char too
    parser = argparse.ArgumentParser(prog_name, prefix_chars='-+')
    parser.add_argument("input", nargs='*', help="the script to be compiled and executed (followed by any command-line arguments")
    parser.add_argument("--preprocess_only",
                        help="only preprocess the input script and not execute it",
                        action="store_true")
    parser.add_argument("--output_preprocessed",
                        help=" output the preprocessed script",
                        action="store_true")
    parser.add_argument("--interactive",
                        help="Executes the script using an interactive internal shell session (experimental)",
                        action="store_true")
    parser.add_argument("-c", "--command",
                        help="Evaluate the following as a script, rather than a file",
                        default=None)
    ## This is not the correct way to parse these, because more than one option can be given together, e.g., -ae
    parser.add_argument("-a",
                        help="Enabling the `allexport` shell option",
                        action="store_true",
                        default=False)
    parser.add_argument("+a",
                        help="Disabling the `allexport` shell option",
                        action="store_false",
                        default=False)    
    ## These two are here for compatibility with respect to bash
    parser.add_argument("-v",
                        help="(experimental) prints shell input lines as they are read",
                        action="store_true")
    parser.add_argument("-x",
                        help="(experimental) prints commands and their arguments as they execute",
                        action="store_true")
    
    config.add_common_arguments(parser)
    args = parser.parse_args()
    config.pash_args = args

    ## Initialize the log file
    config.init_log_file()
    if not config.config:
        config.load_config(args.config_path)

    ## Initialize the graphviz directory
    maybe_init_graphviz_dir(args)

    ## Print all the arguments before they are modified below
    log("Arguments:")
    for arg_name, arg_val in vars(args).items():
        log(arg_name, arg_val)
    log("-" * 40)

    ## TODO: We might need to have a better default (like $0 of pa.sh)
    shell_name = "pash"

    if args.command is not None:
        _, fname = ptempfile()
        with open(fname, 'w') as f:
            f.write(args.command)
        ## If the shell is invoked with -c and arguments after it, then these arguments
        ## need to be assigned to $0, $1, $2, ... and not $1, $2, $3, ...
        if(len(args.input) > 0):
            ## Assign $0
            shell_name = args.input[0]
            args.input = args.input[1:]
        args.input = [fname] + args.input
    elif (len(args.input) > 0):
        shell_name = args.input[0]


    return args, shell_name

def preprocess(ast_objects, config):
    ## This is ids for the temporary files that we will save the IRs in
    irFileGen = FileIdGen()

    ## Preprocess ASTs by replacing AST regions with calls to PaSh's runtime.
    ## Then the runtime will do the compilation and optimization with additional
    ## information.
    preprocessed_asts = replace_ast_regions(ast_objects, irFileGen, config)

    return preprocessed_asts

def shell_env(shell_name: str):
    new_env = os.environ.copy()
    new_env["PASH_TMP_PREFIX"] = config.PASH_TMP_PREFIX
    new_env["pash_shell_name"] = shell_name
    return new_env

## The following two functions need to correspond completely
def bash_prefix_args():
    subprocess_args = ["/usr/bin/env", "bash"]
    ## Add shell specific arguments
    if config.pash_args.a:
        subprocess_args.append("-a")
    else:
        subprocess_args.append("+a")
    if config.pash_args.v:
        subprocess_args.append("-v")
    if config.pash_args.x:
        subprocess_args.append("-x")
    return subprocess_args

def bash_exec_string(shell_name):
    flags = []
    if config.pash_args.a:
        flags.append('-a')
    if config.pash_args.v:
        flags.append('-v')
    if config.pash_args.x:
        flags.append('-x')
    return "exec -a{} bash {} -s $@\n".format(shell_name, " ".join(flags))

def execute_script(compiled_script_filename, command, arguments, shell_name):
    new_env = shell_env(shell_name)
    subprocess_args = bash_prefix_args()
    subprocess_args += ["-c", 'source {}'.format(compiled_script_filename), shell_name] + arguments
    # subprocess_args = ["/usr/bin/env", "bash", compiled_script_filename] + arguments
    log("Executing:", "PASH_TMP_PREFIX={} pash_shell_name={} {}".format(config.PASH_TMP_PREFIX, 
                                                                        shell_name,
                                                                        " ".join(subprocess_args)))
    exec_obj = subprocess.run(subprocess_args, env=new_env, close_fds=False)
    return exec_obj.returncode

if __name__ == "__main__":
    main()
     
