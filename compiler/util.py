from datetime import timedelta
import functools
from typing import Optional, TypeVar, Union, List, Any
TType = TypeVar("TType")
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

def print_time_delta(prefix, start_time, end_time):
    ## Always output time in the log.
    time_difference = (end_time - start_time) / timedelta(milliseconds=1)
    ## If output_time flag is set, log the time
    if (config.OUTPUT_TIME):
        log("{} time:".format(prefix), time_difference, " ms", level=0)
    else:
        log("{} time:".format(prefix), time_difference, " ms")


## This function decorates a function to add the logging prefix (without missing the old prefix)
def logging_prefix(logging_prefix):
    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            old_prefix = config.LOGGING_PREFIX
            config.LOGGING_PREFIX = logging_prefix
            result = func(*args, **kwargs)
            config.LOGGING_PREFIX = old_prefix
            return result
        return wrapper
    return decorator

## This is a wrapper for prints
##
## TODO: Extend the configuration to allow for custom file to output PaSh log. This would
##       allow us to not pollute the .time files.
def log(*args, end='\n', level=1):
    ## If the debug logging level is at least
    ## as high as this log message.
    if (config.DEBUG_LEVEL >= level):
        if(config.LOG_FILE == ""):
            print(config.LOGGING_PREFIX, *args, file=sys.stderr, end=end, flush=True)
        else:
            with open(config.LOG_FILE, "a") as f:
                print(config.LOGGING_PREFIX, *args, file=f, end=end, flush=True)

def ptempfile():
    fd, name = tempfile.mkstemp(dir=config.PASH_TMP_PREFIX)
    ## TODO: Get a name without opening the fd too if possible
    os.close(fd)
    return name

def return_empty_list_if_none_else_itself(arg: Optional[TType]) -> Union[TType, List[Any]]: #list always empty
    if arg is None:
        return []
    else:
        return arg

def return_default_if_none_else_itself(arg: Optional[TType], default: TType) -> TType:
    if arg is None:
        return default
    else:
        return arg

