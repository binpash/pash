import random
import string
from datetime import timedelta
import sys

import config
import tempfile

def flatten_list(lst):
    return [item for sublist in lst for item in sublist]

def pad(lst, index):
    if(index >= len(lst)):
        lst += [None] * (index + 1 - len(lst))
    return lst

def print_time_delta(prefix, start_time, end_time, args):
    if(args.output_time):
        time_difference = (end_time - start_time) / timedelta(milliseconds=1)
        log("{} time:".format(prefix), time_difference, " ms")

## This is a wrapper for prints
##
## TODO: Extend the configuration to allow for custom file to output PaSh log. This would
##       allow us to not pollute the .time files.
def log(*args, end='\n', level=1):
    ## If the debug logging level is at least
    ## as high as this log message.
    if (config.pash_args.debug >= level):
        if(config.pash_args.log_file == ""):
            print(*args, file=sys.stderr, end=end, flush=True)
        else:
            with open(config.pash_args.log_file, "a") as f:
                print(*args, file=f, end=end, flush=True)

def ptempfile():
    return tempfile.mkstemp(dir=config.PASH_TMP_PREFIX)
