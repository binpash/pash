from datetime import timedelta
import os
import sys
import config
import tempfile

def flatten_list(lst):
    return [item for sublist in lst for item in sublist]

def unzip(lst):
    res = [[ i for i, j in lst ],
           [ j for i, j in lst ]]
    return res

def pad(lst, index):
    if(index >= len(lst)):
        lst += [None] * (index + 1 - len(lst))
    return lst

def print_time_delta(prefix, start_time, end_time, args=None):
    ## Always output time in the log.
    time_difference = (end_time - start_time) / timedelta(milliseconds=1)
    ## If output_time flag is set, log the time
    if (config.pash_args.output_time == 1):
        log("{} time:".format(prefix), time_difference, " ms", level=0)
    else:
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
            print(config.LOGGING_PREFIX, *args, file=sys.stderr, end=end, flush=True)
        else:
            with open(config.pash_args.log_file, "a") as f:
                print(config.LOGGING_PREFIX, *args, file=f, end=end, flush=True)

def ptempfile():
    fd, name = tempfile.mkstemp(dir=config.PASH_TMP_PREFIX)
    ## TODO: Get a name without opening the fd too if possible
    os.close(fd)
    return name
