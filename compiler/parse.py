import json
import tempfile
import os
import re
import config
import subprocess
import sys

from util import *

def parse_shell(input_script_path):
    if(not os.path.isfile(input_script_path)):
        log("Error! File:", input_script_path, "does not exist.", level=0)
        exit(1)
    parser_output = subprocess.run([config.PARSER_BINARY, input_script_path], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    if (not parser_output.returncode == 0):
        log(parser_output.stderr)
        parser_output.check_returncode()
    return parser_output.stdout

def from_ir_to_shell(ir_filename):
    printer_output = subprocess.run([config.PRINTER_BINARY, ir_filename], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    printer_output.check_returncode()
    preprocessed_script = printer_output.stdout
    return preprocessed_script

def from_ir_to_shell_file(ir_filename, new_shell_filename):
    preprocessed_script = from_ir_to_shell(ir_filename)
    with open(new_shell_filename, 'w') as new_shell_file:
        new_shell_file.write(preprocessed_script)

