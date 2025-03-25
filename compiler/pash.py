import sys
import os
import subprocess

from ir import *
from pash_graphviz import maybe_init_graphviz_dir
from preprocessor.preprocessor import preprocess
from speculative import util_spec
from util import *
import config
from cli import RunnerParser

LOGGING_PREFIX = "PaSh: "


@logging_prefix(LOGGING_PREFIX)
def main():
    ## Parse arguments
    args, shell_name = parse_args()
    ## If it is interactive we need a different execution mode
    ##
    ## The user can also ask for an interactive mode irregardless of whether pash was invoked in interactive mode.
    if len(args.input) == 0 or args.interactive:
        log("ERROR: --interactive option is not supported!", level=0)
        assert False
    else:
        input_script_path = args.input[0]
        input_script_arguments = args.input[1:]

        ## Preprocess and execute the parsed ASTs
        return_code = preprocess_and_execute_asts(
            input_script_path, args, input_script_arguments, shell_name
        )

        log("-" * 40)  # log end marker
        ## Return the exit code of the executed script
        sys.exit(return_code)


def preprocess_and_execute_asts(
    input_script_path, args, input_script_arguments, shell_name
):
    preprocessed_shell_script = preprocess(input_script_path, args)
    if args.output_preprocessed:
        log("Preprocessed script:")
        log(preprocessed_shell_script)

    ## Write the new shell script to a file to execute
    fname = ptempfile()
    log("Preprocessed script stored in:", fname)
    with open(fname, "w") as new_shell_file:
        new_shell_file.write(preprocessed_shell_script)

    ## 4. Execute the preprocessed version of the input script
    if not args.preprocess_only:
        return_code = execute_script(
            fname, args.command, input_script_arguments, shell_name
        )
    else:
        return_code = 0

    return return_code


def parse_args():
    prog_name = sys.argv[0]
    if "PASH_FROM_SH" in os.environ:
        prog_name = os.environ["PASH_FROM_SH"]
    ## We need to set `+` as a prefix char too
    parser = RunnerParser(prog_name, prefix_chars="-+")
    args = parser.parse_args()
    config.set_config_globals_from_pash_args(args)

    ## Modify the preprocess mode and the partial order file if we are in speculative mode
    if args.speculative:
        log("PaSh is running in speculative mode...")
        args.__dict__["preprocess_mode"] = "spec"
        args.__dict__["partial_order_file"] = util_spec.partial_order_file_path()
        log(" -- Its partial order file will be stored in:", args.partial_order_file)

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
        fname = ptempfile()
        with open(fname, "w") as f:
            f.write(args.command)
        ## If the shell is invoked with -c and arguments after it, then these arguments
        ## need to be assigned to $0, $1, $2, ... and not $1, $2, $3, ...
        if len(args.input) > 0:
            ## Assign $0
            shell_name = args.input[0]
            args.input = args.input[1:]
        args.input = [fname] + args.input
    elif len(args.input) > 0:
        shell_name = args.input[0]

    return args, shell_name


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
        flags.append("-a")
    if config.pash_args.v:
        flags.append("-v")
    if config.pash_args.x:
        flags.append("-x")
    return "exec -a{} bash {} -s $@\n".format(shell_name, " ".join(flags))


def execute_script(compiled_script_filename, command, arguments, shell_name):
    new_env = shell_env(shell_name)
    subprocess_args = bash_prefix_args()
    subprocess_args += [
        "-c",
        "source {}".format(compiled_script_filename),
        shell_name,
    ] + arguments
    # subprocess_args = ["/usr/bin/env", "bash", compiled_script_filename] + arguments
    log(
        "Executing:",
        "PASH_TMP_PREFIX={} pash_shell_name={} {}".format(
            config.PASH_TMP_PREFIX, shell_name, " ".join(subprocess_args)
        ),
    )
    exec_obj = subprocess.run(subprocess_args, env=new_env, close_fds=False)
    return exec_obj.returncode


if __name__ == "__main__":
    main()
