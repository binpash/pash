import sys
import os
import subprocess

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

        ## Preprocess
        preprocess_asts(
            input_script_path, args, input_script_arguments, shell_name
        )

        log("-" * 40)  # log end marker
        sys.exit()


def preprocess_asts(
    input_script_path, args, input_script_arguments, shell_name
):
    preprocessed_shell_script = preprocess(input_script_path, args)
    if args.output_preprocessed:
        log("Preprocessed script:")
        log(preprocessed_shell_script)

    ## Write the new shell script to a file to execute
    ## Use --output path if provided, otherwise create temp file
    if args.output:
        fname = args.output
    else:
        fname = ptempfile()

    log("Preprocessed script stored in:", fname)
    with open(fname, "wb") as new_shell_file:
        preprocessed_shell_script = preprocessed_shell_script.encode(
            "utf-8", errors="replace"
        )
        new_shell_file.write(preprocessed_shell_script)

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




if __name__ == "__main__":
    main()
