import sys
import os
import subprocess
import argparse
import json
import logging

from preprocessor import preprocess
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
            help="the script file to be preprocessed",
        )
        self.add_argument(
            "--output",
            help="path where the preprocessed script will be saved",
            required=True,
        )
        self.add_argument(
            "--bash",
            help="(experimental) interpret the input as a bash script file",
            action="store_true",
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
    args = parse_args()
    input_script_path = args.input

    ## Preprocess
    preprocess_asts(
        input_script_path, args
    )

    log("-" * 40)  # log end marker

def preprocess_asts(
    input_script_path, args
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
    parser = Parser(prog_name)
    args = parser.parse_args()
    config_from_args(args)

    ## Modify the preprocess mode and the partial order file if we are in speculative mode
    if args.speculative:
        log("PaSh is running in speculative mode...")
        args.__dict__["preprocess_mode"] = "spec"
        args.__dict__["partial_order_file"] = util_spec.partial_order_file_path()
        log(" -- Its partial order file will be stored in:", args.partial_order_file)

    ## Print all the arguments
    log("Arguments:")
    for arg_name, arg_val in vars(args).items():
        log(arg_name, arg_val)
    log("-" * 40)

    return args




if __name__ == "__main__":
    main()
