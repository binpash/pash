import sys
import os
import re

unix50_dir = sys.argv[1]
intermediaries_dir = sys.argv[2]
input_size_increase = int(sys.argv[3])

## TODO: Make this variable depend on the input size increase
MAXIMUM_INPUT_SIZE = 10000000000 # 10 GB

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
input_file_name_assignments = {line.split("=")[0] : line.split("=")[1] for line in input_lines}
input_file_names = input_file_name_assignments.values()

for input_file in input_file_names:
    input_file_path = os.path.join(unix50_dir, input_file)
    generated_input_file_path = os.path.join(generated_inputs_dir, input_file)
    with open(input_file_path) as file:
        input_file_data = file.read()

    # print("Input:", input_file_path, "size:", len(input_file_data))
    # print("Increased size:", input_size_increase * len(input_file_data))
    num_iterations = min(input_size_increase, MAXIMUM_INPUT_SIZE // len(input_file_data))
    with open(generated_input_file_path, "w") as file:
        for _ in range(num_iterations):
            file.write(input_file_data)

## Extract each pipeline of unix50
unix50_pipelines = [line for line in unix50_script_lines if line.startswith("cat ")]
p = re.compile('ab*')

for i, unix50_pipeline in enumerate(unix50_pipelines):

    ## Find the input's name
    input_name = re.search('IN[0-9]+', unix50_pipeline).group(0)
    unix50_normalized_pipeline = re.sub('IN[0-9]+', 'IN', unix50_pipeline)

    ## Generate the script
    pipeline_file_path = os.path.join(intermediaries_dir, "unix50_pipeline_{0:03d}.sh".format(i))
    with open(pipeline_file_path, "w") as file:
        file.write(unix50_normalized_pipeline + "\n")

    ## And its environment
    pipeline_env_file_path = os.path.join(intermediaries_dir, "unix50_pipeline_{0:03d}_env.sh".format(i))

    input_file_assignment_path = os.path.join(generated_inputs_dir, input_file_name_assignments[input_name])
    environment_data = "IN={}\n".format(input_file_assignment_path)
    with open(pipeline_env_file_path, "w") as file:
        file.write(environment_data)
