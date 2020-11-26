import random
import string
from datetime import timedelta
import sys

def flatten_list(lst):
    return [item for sublist in lst for item in sublist]

def print_time_delta(prefix, start_time, end_time, args):
    if(args.output_time):
        time_difference = (end_time - start_time) / timedelta(milliseconds=1)
        print("{} time:".format(prefix), time_difference, " ms", file=sys.stderr)

def get_random_string(length=8):
    letters = string.ascii_lowercase
    result_str = ''.join(random.choice(letters) for i in range(length))
    return result_str