import os
import shlex
import subprocess
import sys

GIT_TOP_CMD = [ 'git', 'rev-parse', '--show-toplevel', '--show-superproject-working-tree']
if 'PASH_TOP' in os.environ:
    PASH_TOP = os.environ['PASH_TOP']
else:
    PASH_TOP = subprocess.run(GIT_TOP_CMD, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True).stdout.rstrip()

sys.path.append(os.path.join(PASH_TOP, "compiler"))

## First set the tmp prefix because config needs it
os.environ['PASH_TMP_PREFIX'] = '/tmp'

import config


def find_bug_with_first_output():
    with open("output.out") as f:
        data = f.read()

        lines = data.split("\n")

        # 1031
        # 1315
        # 1397
        ## 1420
        target_data = "\n".join(lines[1431:1432])
        # target_data = data
        print("data:", target_data)
        tokens = shlex.split(target_data, comments=True, posix=True)
        print("Done")
        print(tokens)

def find_bug_with_second_output():
    with open("output9.out") as f:
        data = f.read()

        lines = data.split("\n")


        target_data = "\n".join(lines[3453:3454])
        # target_data = data
        print("data:", target_data)
        tokens = shlex.split(target_data, comments=True, posix=False)
        print("Done")
        print(tokens)

# find_bug_with_second_output()

def test_var_file_read():
    ## Set some config state
    class Object:
        pass

    config.pash_args = Object()
    config.pash_args.debug = 1
    config.pash_args.log_file = ""
    filename = "test-shlex-aux.sh"
    
    config.read_vars_file(filename)


test_var_file_read()