import random
import string
from datetime import timedelta
import sys

import config

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
def log(*args, end='\n'):
    if(config.pash_args.log_file == ""):
        if (not config.pash_args.debug == "0"):
            print(*args, file=sys.stderr, end=end)
    else:
        with open(config.pash_args.log_file, "a") as f:
            print(*args, file=f, end=end)
    

def get_random_string(length=8):
    letters = string.ascii_lowercase
    result_str = ''.join(random.choice(letters) for i in range(length))
    return result_str
