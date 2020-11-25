import json
import tempfile
import os
import re
import config
import subprocess


def parse_shell(input_script_path):
    parser_output = subprocess.run([config.PARSER_BINARY, input_script_path], capture_output=True, text=True)
    if (not parser_output.returncode == 0):
        print(parser_output.stderr)
        parser_output.check_returncode()
    return parser_output.stdout

def from_ir_to_shell(ir_filename):
    printer_output = subprocess.run([config.PRINTER_BINARY, ir_filename], capture_output=True, text=True)
    printer_output.check_returncode()
    preprocessed_script = printer_output.stdout
    return preprocessed_script

def from_ir_to_shell_file(ir_filename, new_shell_filename):
    preprocessed_script = from_ir_to_shell(ir_filename)
    with open(new_shell_filename, 'w') as new_shell_file:
        new_shell_file.write(preprocessed_script)


##
## Obsolete code
##

def shell_to_ir(script_contents, new_file=False):
    with tempfile.NamedTemporaryFile(mode="w+", prefix="sh", suffix="temp", delete=False) as file:
    # with open("tmp-DELETE.sh", "w+") as file:
        file.write(script_contents)
        file.close()
    try:
        return shell_file_to_ir(file.name, new_file)
    finally:
        os.remove(file.name)

def ir_to_shell(ir, new_file=False):
    with tempfile.NamedTemporaryFile(mode="w+", prefix="ir", suffix="temp", delete=False) as file:
        file.write("\n".join([json.dumps(ast) for ast in ir]))
        file.close()
    try:
        return ir_file_to_shell(file.name, new_file)
    finally:
        os.remove(file.name)

def shell_file_to_ir(filename, new_file=False):
    parser_output = subprocess.run([config.PARSER_BINARY, filename], capture_output=True, text=True)
    if (not parser_output.returncode == 0):
        print(parser_output.stderr)
        parser_output.check_returncode()
    if (new_file):
        with open(new_file, 'w') as f:
            f.write(parser_output.stdout)
    else:
        s = parser_output.stdout
        return [json.loads(l) for l in s.strip().split("\n")]

def ir_file_to_shell(filename, new_file=False):
    printer_output = subprocess.run([config.PRINTER_BINARY, filename], capture_output=True, text=True)
    if (not printer_output.returncode == 0):
        print(printer_output.stderr)
        printer_output.check_returncode()
    if (new_file):
        with open(new_file, 'w') as f:
            f.write(printer_output.stdout)
    else:
        s = printer_output.stdout
