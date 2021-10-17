import argparse
import traceback

from annotations import *
import config
import pash_runtime
from util import *

##
## A Daemon responding to requests for compilation
##
## Note: Not an actual daemon with the strict Unix sense
##

## TODO: Rename the pash_runtime to pash_compiler and this to pash_daemon

## TODO: Should we maybe use sockets instead of fifos?

## TODO: Fix the daemon logging.

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="the input fifo from which the daemon will read its input")
    parser.add_argument("output", help="the output fifo to which the daemon will write its output")
    config.add_common_arguments(parser)
    args, unknown_args = parser.parse_known_args()

    return args

## Initialize the daemon
def init():
    args = parse_args()
    config.pash_args = args

    ## Load the configuration
    if not config.config:
        config.load_config(args.config_path)
    
    ## Load annotations
    config.annotations = load_annotation_files(config.config['distr_planner']['annotations_dir'])

    pash_runtime.runtime_config = config.config['distr_planner']

    return args

def success_response(string):
    return f'OK: {string}\n'

def error_response(string):
    return f'ERROR: {string}\n'

def parse_command(input):
    ## TODO: Improve the way parsing happens plz :')
    if(input.startswith("Compile:")):
        return compile(input)
    else:
        return error_response(f'Unsupported command: {input}')

## TODO: Improve the way parsing happens plz :') At the moment this will fail with : in file etc
def parse_compile_command(input):
    try:
        components = input.rstrip().split("|")
        compiled_script_file = components[0].split(":")[1]
        var_file = components[1].split(":")[1]
        input_ir_file = components[2].split(":")[1]
        return compiled_script_file, var_file, input_ir_file
    except:
        raise Exception(f'Parsing failure for line: {input}')

def compile(input):
    compiled_script_file, var_file, input_ir_file = parse_compile_command(input)
    
    ## Read any shell variables files if present
    config.read_vars_file(var_file)

    ## Call the main procedure
    pash_runtime.compile_optimize_output_script(input_ir_file, compiled_script_file, config.pash_args)

    return success_response(f'{compiled_script_file} {var_file} {input_ir_file}')

class Scheduler:
    def __init__(self):
        self.input_fids = set()
        self.output_fids = set()
        self.process_fids = {} # map process_id -> (input_fids, output_fids)
        self.next_id = 0
        self.waiting_process = (None, None) # (process_id, compilation_response)
        self.running_procs = 0
        self.unsafe_running = False
        self.done = False
        
    def check_fids_safety(self, process_id):
        proc_input_fids, proc_output_fids = self.process_fids[process_id]
        all_proc_fids = proc_input_fids.union(proc_output_fids)
        if self.output_fids.intersection(all_proc_fids) or self.input_fids.intersection(proc_output_fids):
            return False
        return True

    def compile_and_add(self, compiled_script_file, var_file, input_ir_file):
        self.wait_unsafe()
        process_id = self.get_next_id()
        # optimized_ast_or_ir = pash_runtime.compile_ir(input_ir_file, config.pash_args)
        compile_success = False
        run_parallel = False
        ## Read any shell variables files if present
        config.read_vars_file(var_file)
        ast_or_ir = pash_runtime.compile_ir(input_ir_file, compiled_script_file, config.pash_args)

        if ast_or_ir == None:
            response = error_response(f'{process_id} failed to compile')
            self.waiting_process = (process_id, response)
            return self.wait_for_all()
        else:
            response = success_response(f'{process_id} {compiled_script_file} {var_file} {input_ir_file}')
            compile_success = True

            proc_input_fids = set(map(lambda out: str(out.resource) if str(out.resource) != "None" else out, ast_or_ir.all_input_fids()))
            proc_output_fids = set(map(lambda out: str(out.resource) if str(out.resource) != "None" else out, ast_or_ir.all_output_fids()))
        
            self.process_fids[process_id] = (proc_input_fids, proc_output_fids)
            run_parallel = self.check_fids_safety(process_id)
            if run_parallel:
                self.input_fids = self.input_fids.union(proc_input_fids)
                self.output_fids = self.output_fids.union(proc_output_fids)
                self.running_procs += 1
            else:
                self.waiting_process = (process_id, response)
                return self.wait_for_all()

            self.fout.write(response)
            

    def remove_process(self, process_id):
        if process_id in self.process_fids:
            proc_input_fids, proc_output_fids = self.process_fids[process_id]
            self.input_fids -= proc_input_fids
            self.output_fids -= proc_output_fids
            del self.process_fids[process_id]
            
        self.running_procs -= 1
        if self.running_procs == 0:
            self.unsafe_running = False

    def is_idle(self):
        return len(self.fids_in_use) + len(self.process_fids) == 0

    def get_next_id(self):
        self.next_id += 1
        return self.next_id

    def wait_for_all(self):
        while self.running_procs > 0:
            with open(self.in_filename) as fin:
                input_cmd = fin.read()
                # must be exit command
                if (input_cmd.startswith("Exit:")):
                    self.remove_process(int(input_cmd.split(":")[1]))
                else:
                    raise Exception(f"Command should be exit but it was {input_cmd}")
                # print(f"waitall {input_cmd}")
        # finished waiting -> execute waiting process and wait for it to finish
        self.run_waiting_process()

    def run_waiting_process(self):
        waiting_process_id, response = self.waiting_process
        if not waiting_process_id:
            return
        self.fids_in_use = set()
        self.waiting_process = (None, None)
        self.fout.write(response)
        self.running_procs += 1
        self.unsafe_running = True

    def wait_unsafe(self):
        if self.unsafe_running:
            with open(self.in_filename) as fin:
                input_cmd = fin.read()
                if (not input_cmd.startswith("Exit:")):
                    raise Exception(f"Command should be exit but it was {input_cmd}")
                #print(f"unsafe {input_cmd}")
            self.unsafe_running = False
            self.running_procs -= 1

    def parse_and_run_cmd(self, input_cmd):
        if(input_cmd.startswith("Compile")):
            compiled_script_file, var_file, input_ir_file = parse_compile_command(input_cmd)
            self.compile_and_add(compiled_script_file, var_file, input_ir_file)
        elif (input_cmd.startswith("Exit")):
            #print(input_cmd.split(":")[1], file=sys.stderr)
            self.remove_process(int(input_cmd.split(":")[1]))
        elif (input_cmd.startswith("Done")):
            self.wait_for_all()
            self.fout.write("All finished")
            self.done = True
        else:
            self.fout.write(error_response(f'Unsupported command: {input_cmd}'))

    def run(self, in_filename, out_filename):
        while not self.done:
            self.in_filename = in_filename
            self.out_filename = out_filename
            ## Process a single request
            with open(in_filename) as fin:
                input_cmd = fin.read()
                if not input_cmd:
                    continue
                input_cmd = input_cmd.rstrip()
                # print(input_cmd)
                if (not input_cmd.startswith("Exit:")):
                    self.fout = open(out_filename, "w")

                self.parse_and_run_cmd(input_cmd)
                
                if (not input_cmd.startswith("Exit:")):
                    self.fout.flush()
                    self.fout.close()
                # if self.running_procs == 0 and self.wait_for_all:
                #     self.wait_for_all = False
                #     if self.waiting_process[0]:
                #         waiting_process_id, response = self.waiting_process
                #         self.running_procs += 1
                #         self.fids_in_use = set()
                #         # Should always be true for now
                #         self.waiting_process = (None, None)
                #         fout.write(response)
                #         fout.flush()



def main():
    args = init()
    scheduler = Scheduler()
    scheduler.run(in_filename=args.input, out_filename=args.output)

if __name__ == "__main__":
    main()
#TODO
### Flag to turn off, abstaction, interferancem