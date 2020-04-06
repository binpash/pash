import sys
import os
from shutil import copyfile

unix50_dir = sys.argv[1]
intermediaries_dir = sys.argv[2]
input_size_increase = int(sys.argv[3])

## Make the generated inputs dir
generated_inputs_dir = os.path.join(unix50_dir, "inputs")
if not os.path.exists(generated_inputs_dir):
    os.makedirs(generated_inputs_dir)

## Read the unix50 script to find inputs and separate pipelines
unix50_script = os.path.join(unix50_dir, "unix50.sh")
with open(unix50_script) as file:
    unix50_script_data = file.read()

## Generate bigger inputs
unix50_script_lines = unix50_script_data.split("\n")
input_lines = [line for line in unix50_script_lines if line.startswith("IN")]
input_file_names = [line.split("=")[1] for line in input_lines]

for input_file in input_file_names:
    input_file_path = os.path.join(unix50_dir, input_file)
    generated_input_file_path = os.path.join(generated_inputs_dir, input_file)
    with open(input_file_path) as file:
        input_file_data = file.read()

    with open(generated_input_file_path, "w") as file:
        file.write(input_file_data * input_size_increase)


## TODO: Extract each pipeline of unix50
