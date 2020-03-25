import sys
import os
from shutil import copyfile

input_dir = sys.argv[1]
name_of_script = sys.argv[2]
number_of_inputs = int(sys.argv[3])
output_dir = sys.argv[4]

## This script takes a microbenchmark script as input, finds the $IN
## occurence in it and then generates an intermediary script with many
## $INs in its place.

input_script = os.path.join(input_dir, name_of_script + ".sh")
output_script = os.path.join(output_dir, '{}_{}_seq.sh'.format(name_of_script, number_of_inputs))
input_env = os.path.join(input_dir, name_of_script + "_env.sh")
output_env = os.path.join(output_dir, '{}_{}_env.sh'.format(name_of_script, number_of_inputs))
input_funs = os.path.join(input_dir, name_of_script + "_funs.sh")
output_funs = os.path.join(output_dir, '{}_{}_funs.sh'.format(name_of_script, number_of_inputs))

## Generate the sequential script
with open(input_script) as file:
    input_script_data = file.read()

output_script_data = input_script_data.replace(' $IN', ' $IN' * number_of_inputs, 1)

with open(output_script, "w") as file:
    file.write(output_script_data)

## Copy the environment script
copyfile(input_env, output_env)

## Copy the funs file (if it exists)
if os.path.exists(input_funs):
    copyfile(input_funs, output_funs)
