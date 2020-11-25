import sys
import os
import subprocess
from shutil import copyfile, rmtree

def find_and_split_inputs(input_env, output_dir, number_of_inputs):
    ## Generate the environment file (by splitting inputs)
    with open(input_env) as file:
        input_env_data = file.read()

    ## Find the input file name
    input_env_lines = input_env_data.split('\n')
    input_vars = [line.split('=')[1] for line in input_env_lines if line.startswith('IN=')]
    assert(len(input_vars) == 1)
    input_file_name = input_vars[0]
    # print(input_file_name)

    ## Split input files
    new_input_files = split_inputs(output_dir, number_of_inputs, input_file_name)

    return new_input_files

def split_inputs(output_dir, number_of_inputs, input_file_name):
    ## Make a directory to store the split parts (and delete what was previously there)
    split_input_directory = os.path.join(output_dir, 'split_inputs')
    if os.path.exists(split_input_directory):
        rmtree(split_input_directory)

    os.makedirs(split_input_directory)

    ## Split it into parts
    # stream = os.popen('wc -l {}'.format(input_file_name))
    # wc_output = stream.read()
    # input_file_n_lines = wc_output.split()[0]
    # print(input_file_n_lines)
    split_file_prefix = os.path.join(split_input_directory, 'input-chunk-')
    stream = subprocess.run(['gsplit',
                             '-n l/{}'.format(number_of_inputs),
                             '-d',
                             '{}'.format(input_file_name),
                             '{}'.format(split_file_prefix)],
                             check=True)

    new_input_files = [os.path.join(split_input_directory, f) for f in os.listdir(split_input_directory)
                       if os.path.isfile(os.path.join(split_input_directory, f))]

    new_input_files.sort()
    # print(new_input_files)
    return new_input_files

def list_split_inputs(output_dir):
    split_input_directory = os.path.join(output_dir, 'split_inputs')
    assert(os.path.exists(split_input_directory))

    new_input_files = [os.path.join(split_input_directory, f) for f in os.listdir(split_input_directory)
                       if os.path.isfile(os.path.join(split_input_directory, f))]

    new_input_files.sort()
    # print(new_input_files)
    return new_input_files


def generate_env_file(input_env, output_env, new_input_files):
    ## Generate the environment file (by splitting inputs)
    with open(input_env) as file:
        input_env_data = file.read()

    ## Find the input file name
    input_env_lines = input_env_data.split('\n')

    ## Save the new environment file accordingly
    no_input_env_lines = [line for line in input_env_lines if not line.startswith('IN=')]
    new_input_vars = ['IN{}={}'.format(i, in_file_name)
                      for i, in_file_name in enumerate(new_input_files)]
    output_env_data = "\n".join(new_input_vars + no_input_env_lines)
    with open(output_env, "w") as file:
        file.write(output_env_data)

## Replace $IN with all different $INis (one for each input file)
def replace_in_variable(data, new_input_files):
    new_data = data.replace(' $IN', ' ' + ' '.join(['$IN{}'.format(i)
                                                    for i in range(len(new_input_files))]))
    return new_data

def generate_seq_script(input_script, output_script, new_input_files):
    ## Generate the sequential script
    with open(input_script) as file:
        input_script_data = file.read()

    output_script_data = replace_in_variable(input_script_data, new_input_files)

    with open(output_script, "w") as file:
        file.write(output_script_data)

def main():
    input_dir = sys.argv[1]
    name_of_script = sys.argv[2]
    number_of_inputs = int(sys.argv[3])
    output_dir = sys.argv[4]

    try:
        env_suffix = sys.argv[5]
    except:
        env_suffix = "env"

    ## This script takes a microbenchmark script as input, finds the $IN
    ## occurence in it and then generates an intermediary script with many
    ## $INs in its place.
    input_script = os.path.join(input_dir, name_of_script + ".sh")
    output_script = os.path.join(output_dir, '{}_{}_seq.sh'.format(name_of_script, number_of_inputs))
    input_env = os.path.join(input_dir, name_of_script + "_{}.sh".format(env_suffix))
    output_env = os.path.join(output_dir, '{}_{}_env.sh'.format(name_of_script, number_of_inputs))
    input_funs = os.path.join(input_dir, name_of_script + "_funs.sh")
    output_funs = os.path.join(output_dir, '{}_{}_funs.sh'.format(name_of_script, number_of_inputs))

    ## Find and split input files given the environment file
    new_input_files = find_and_split_inputs(input_env, output_dir, number_of_inputs)

    ## Generate new environment file
    generate_env_file(input_env, output_env, new_input_files)

    ## Copy the funs file (if it exists)
    if os.path.exists(input_funs):
        copyfile(input_funs, output_funs)

    generate_seq_script(input_script, output_script, new_input_files)


if __name__ == "__main__":
    main()

