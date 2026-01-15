import sys
import os
import subprocess
import argparse
import json
import logging

from preprocessor.preprocessor import preprocess
from speculative import util_spec
from util import *

LOGGING_PREFIX = "PaSh-Preprocessor: "

def config_from_args(pash_args):
    ## Also set logging here
    # Format logging
    # ref: https://docs.python.org/3/library/logging.html#formatter-objects
    ## TODO: When we add more logging levels bring back the levelname+time
    if pash_args.log_file == "":
        logging.basicConfig(format="%(message)s")
    else:
        logging.basicConfig(
            format="%(message)s",
            filename=f"{os.path.abspath(pash_args.log_file)}",
            filemode="w",
        )

    # Set debug level
    if pash_args.debug == 1:
        logging.getLogger().setLevel(logging.INFO)
    elif pash_args.debug >= 2:
        logging.getLogger().setLevel(logging.DEBUG)

class Parser(argparse.ArgumentParser):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.add_argument(
            "-d",
            "--debug",
            type=int,
            help="configure debug level; defaults to 0",
            default=0,
        )
        self.add_argument(
            "--log_file",
            help="configure where to write the log; defaults to stderr.",
            default="",
        )
        self.add_argument(
            "input",
            nargs="*",
            help="the script to be compiled and executed (followed by any command-line arguments",
        )
        self.add_argument(
            "--output",
            help="path where the preprocessed script will be saved (used with preprocessing mode)",
            required=True,
        )
        self.add_argument(
            "--bash",
            help="(experimental) interpret the input as a bash script file",
            action="store_true",
        )
        self.add_argument(
            "-c",
            "--command",
            help="evaluate the following as a script, rather than a file",
            default=None,
        )
        self.add_argument(
            "--speculative",
            help="(experimental) use the speculative execution preprocessing and runtime",
            action="store_true",
            default=False,
        )

        self.set_defaults(preprocess_mode="pash")

@logging_prefix(LOGGING_PREFIX)
def main():
    ## Parse arguments
    args, shell_name = parse_args()    
    input_script_path = args.input[0]

    ## Preprocess
    preprocess_asts(
        input_script_path, args, shell_name
    )

    log("-" * 40)  # log end marker

def preprocess_asts(
    input_script_path, args, shell_name
):
    preprocessed_shell_script = preprocess(input_script_path, args)

    ## Write the new shell script to a file to execute
    fname = args.output

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
    parser = Parser(prog_name, prefix_chars="-+")
    ## Use parse_known_args to allow pa.sh to pass arguments meant for other components
    ## (e.g., compilation server arguments like --width, --no_optimize)
    args, unknown_args = parser.parse_known_args()
    config_from_args(args)
    if unknown_args:
        log(f"Ignoring unknown arguments: {' '.join(unknown_args)}", level=1)

    ## Modify the preprocess mode and the partial order file if we are in speculative mode
    if args.speculative:
        log("PaSh is running in speculative mode...")
        args.__dict__["preprocess_mode"] = "spec"
        args.__dict__["partial_order_file"] = util_spec.partial_order_file_path()
        log(" -- Its partial order file will be stored in:", args.partial_order_file)

    ## Print all the arguments before they are modified below
    log("Arguments:")
    for arg_name, arg_val in vars(args).items():
        log(arg_name, arg_val)
    log("-" * 40)

    ## TODO: Can the following code be removed?

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
