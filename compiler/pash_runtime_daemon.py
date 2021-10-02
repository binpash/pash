import argparse

import config
# import pash_runtime

##
## A Daemon responding to requests for compilation
##
## Note: Not an actual daemon with the strict Unix sense
##

## TODO: Rename the pash_runtime to pash_compiler and this to pash_daemon

## TODO: Should we maybe use sockets instead of fifos?

## TODO: Fix the daemon logging.

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="the input fifo from which the daemon will read its input")
    parser.add_argument("output", help="the output fifo to which the daemon will write its output")
    config.add_common_arguments(parser)
    args, unknown_args = parser.parse_known_args()

    ## Print all the arguments before they are modified below
    print("Arguments:")
    for arg_name, arg_val in vars(args).items():
        print(arg_name, arg_val)
    print("-" * 40)
    print("Unknown rest arguments:")
    print(unknown_args)
    print("-" * 40)

    return args

## Initialize the daemon
def init():
    args = parse_args()
    return args

## TODO: Implement actual commands here
def parse_line(line):
    return line

def main():
    args = init()

    while True:
        ## Process a single request
        with open(args.input) as fin, open(args.output, "w") as fout:
            for line in fin.readlines():
                ret = parse_line(line)
                fout.write("OK: " + ret)

if __name__ == "__main__":
    main()