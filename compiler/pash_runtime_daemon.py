import argparse
import pexpect
import signal
import subprocess
import sys
import traceback

from annotations import *
import config
import pash_runtime
from util import *

##
# A Daemon responding to requests for compilation
##
# Note: Not an actual daemon with the strict Unix sense
##

# TODO: Rename the pash_runtime to pash_compiler and this to pash_daemon

# TODO: Should we maybe use sockets instead of fifos?

# TODO: Fix the daemon logging.


def handler(signum, frame):
    log("Signal:", signum, "caught")
    shutdown()

signal.signal(signal.SIGTERM, handler)

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "input", help="the input fifo from which the daemon will read its input")
    parser.add_argument(
        "output", help="the output fifo to which the daemon will write its output")
    config.add_common_arguments(parser)
    args, unknown_args = parser.parse_known_args()

    return args

# Initialize the daemon


def init():
    ## Set the logging prefix
    config.LOGGING_PREFIX = "Daemon: "
    
    args = parse_args()
    config.pash_args = args

    # Load the configuration
    if not config.config:
        config.load_config(args.config_path)

    # Load annotations
    config.annotations = load_annotation_files(
        config.config['distr_planner']['annotations_dir'])

    pash_runtime.runtime_config = config.config['distr_planner']

    if config.pash_args.expand_using_bash_mirror:
        ## Initialize a bash that is used for expanding
        ##
        ## TODO: Alternatively, we could set up a communication with the original bash 
        ## (though this makes it difficult to have concurrent compilations and execution)
        ## TODO: We actually need to discuss which arch is better.
        bash_mirror = init_bash_mirror_subprocess()

        ## Is it OK to save it in config?
        config.bash_mirror = bash_mirror

    return args

def init_bash_mirror_subprocess():
    ## Spawn a bash process to ask it for expansions
    p = pexpect.spawn('/usr/bin/env', ['bash', '-i'], 
                      encoding='utf-8',
                      echo=False)
    ## If we are in debug mode also log the bash's output
    if (config.pash_args.debug >= 1):
        _, file_to_save_output = ptempfile()
        log("bash mirror log saved in:", file_to_save_output)
        fout = open(file_to_save_output, "w")
        p.logfile = fout
    return p

def success_response(string):
    return f'OK: {string}\n'


def error_response(string):
    return f'ERROR: {string}\n'


class Scheduler:
    """ Takes care of running processes in parallel if there is no conflict. 
    The scheduler relies on the fact that process will wait for a compilation response.
    This allows it to control wether to allow the next process to run or wait for all other process.
    Flow:
        input cmd -> 
                    |   Compile -> 
                            1- Try compiling the pipeline
                            2- Wait for any unsafe processes to finish
                            3- Check compilation for success and any conficts 
                                - no side effects -> allow to run in parallel by sending a response
                                - failed or conflict -> wait for all process to exit then run this process in unsafe mode

                    |   Exit process_id -> remove process_id from the list of running processes

                    |   Done -> no more pipelines -> wait for all processes to finish and exit

    Notes:
        The current design relies on the data being written atomicly to the pipe for the size of our input. This allows as to use only one pipe instead of a pipe per process.
        Source: https://www.gnu.org/software/libc/manual/html_node/Pipe-Atomicity.html
    """

    def __init__(self):
        self.input_resources = set()
        self.output_resources = set()
        self.process_resources = {}  # map process_id -> (input_resources, output_resources)
        self.next_id = 0
        self.running_procs = 0
        self.unsafe_running = False
        self.done = False
        self.cmd_buffer = ""
        self.reader = None
        self.reader_pipes_are_blocking = True

    def check_resources_safety(self, process_id):
        proc_input_resources, proc_output_resources = self.process_resources[process_id]
        all_proc_resources = proc_input_resources.union(proc_output_resources)
        if self.output_resources.intersection(all_proc_resources) or self.input_resources.intersection(proc_output_resources):
            return False
        return True

    def compile_and_add(self, compiled_script_file, var_file, input_ir_file):
        process_id = self.get_next_id()
        run_parallel = False
        compile_success = False

        # Read any shell variables files if present
        config.read_vars_file(var_file)
        if config.pash_args.expand_using_bash_mirror:
            ## Update the bash mirror with the new variables
            config.update_bash_mirror_vars(var_file)
            ## TODO: Maybe we also need to update current directory of bash mirror for file-based expansion?

            ## Clean up the variable cache
            config.reset_variable_cache()

        ast_or_ir = pash_runtime.compile_ir(
            input_ir_file, compiled_script_file, config.pash_args)
        
        self.wait_unsafe()
        if ast_or_ir != None:
            compile_success = True

            proc_input_resources = set(map(lambda out: str(out.resource) if str(
                out.resource) != "None" else out, ast_or_ir.all_input_fids()))
            proc_output_resources = set(map(lambda out: str(out.resource) if str(
                out.resource) != "None" else out, ast_or_ir.all_output_fids()))

            self.process_resources[process_id] = (proc_input_resources, proc_output_resources)

            run_parallel = self.check_resources_safety(process_id)
            if run_parallel:
                self.input_resources = self.input_resources.union(proc_input_resources)
                self.output_resources = self.output_resources.union(proc_output_resources)

        
        if not run_parallel:
            self.wait_for_all()
            
        if compile_success:
            response = success_response(
                f'{process_id} {compiled_script_file} {var_file} {input_ir_file}')
        else:
            response = error_response(f'{process_id} failed to compile')
            self.unsafe_running = True

        self.running_procs += 1
        return response

    def remove_process(self, process_id):
        log("The following process exited:", process_id)
        if process_id in self.process_resources:
            del self.process_resources[process_id]
            # TODO: Should be improved to not rebuild inputs and outputs from scratch maybe use counters
            self.input_resources = set().union(*[input_resources for input_resources, _ in self.process_resources.values()])
            self.output_resources = set().union(*[output_resources for _, output_resources in self.process_resources.values()])

        self.running_procs -= 1
        if self.running_procs == 0:
            self.unsafe_running = False

    def get_next_id(self):
        self.next_id += 1
        return self.next_id

    def wait_for_all(self):
        log("Waiting for all process to finish:", self.running_procs)
        while self.running_procs > 0:
            input_cmd = self.get_input()
            # must be exit command or something is wrong
            if (input_cmd.startswith("Exit:")):
                self.remove_process(int(input_cmd.split(":")[1]))
            else:
                raise Exception(
                    f"Command should be exit but it was {input_cmd}")
        self.unsafe_running = False

    def wait_unsafe(self):
        log("Unsafe running:", self.unsafe_running)
        if self.unsafe_running:
            assert(self.running_procs == 1)
            self.wait_for_all()
            self.unsafe_running = False

    def parse_and_run_cmd(self, input_cmd):
        if(input_cmd.startswith("Compile")):
            compiled_script_file, var_file, input_ir_file = self.__parse_compile_command(
                input_cmd)
            response = self.compile_and_add(compiled_script_file, var_file, input_ir_file)
            ## Send output to the specific command
            self.send_output(response)
        elif (input_cmd.startswith("Exit")):
            self.remove_process(int(input_cmd.split(":")[1]))
        elif (input_cmd.startswith("Done")):
            self.wait_for_all()
            ## We send output to the top level pash process
            ## to signify that we are done
            self.send_output("All finished")
            self.done = True
        else:
            self.send_output(error_response(
                f'Unsupported command: {input_cmd}'))

    ## This method calls the reader to get an input
    def get_input(self):
        cmd = ""
        if not self.reader_pipes_are_blocking:
            while not cmd:
                cmd = self.reader.get_next_cmd()
        else:
            cmd = self.reader.get_next_cmd()
        return cmd

    ## This method sends output through the output fifo of the daemon
    def send_output(self, message):
        self.fout = open(self.out_filename, "w")
        self.fout.write(message)
        self.fout.flush()
        self.fout.close()

    def __parse_compile_command(self, input):
        try:
            components = input.rstrip().split("|")
            compiled_script_file = components[0].split(":")[1]
            var_file = components[1].split(":")[1]
            input_ir_file = components[2].split(":")[1]
            return compiled_script_file, var_file, input_ir_file
        except:
            raise Exception(f'Parsing failure for line: {input}')

    def run(self, in_filename, out_filename):
        self.reader = UnixPipeReader(in_filename, self.reader_pipes_are_blocking)
        self.out_filename = out_filename
        while not self.done:
            # Process a single request
            input_cmd = self.get_input()

            ## Parse the command (potentially also sending a response)
            self.parse_and_run_cmd(input_cmd)
        
        self.reader.close()
        shutdown()


class UnixPipeReader:
    def __init__(self, in_filename, blocking = True):
        self.in_filename = in_filename
        self.buffer = ""
        self.blocking = blocking
        if not self.blocking:
            # Non blocking mode shouldn't be used in production. It's only used experimentally.
            log("Reader initialized in non-blocking mode")
            self.fin = open(self.in_filename)
        else:
            log("Reader initialized in blocking mode")

    def get_next_cmd(self):
        """
        This method return depends on the reading mode. In blocking mode this method will
        return the next full command and if there is no command it will wait until a full command is recieved.
        In non blocking mode it would either a full command or an empty string if a full command isn't available yet.
        This command keeps a state of the remaining data which is used in each subsequent call to this method.
        """
        input_buffer = ""
        if self.buffer:
            # Don't wait on fin if cmd buffer isn't empty
            log("Reader buffer isn't empty. Using it instead of reading new data for the next command")
            input_buffer = self.buffer
        else:
            log("Reader buffer is empty. Reading new data from input fifo")
            if self.blocking:
                with open(self.in_filename) as fin:
                    # This seems to be necessary for reading the full data. 
                    # It seems like slower/smaller machines might not read the full data in one read
                    while True:
                        data = fin.read()
                        if len(data) == 0:
                            break
                        input_buffer += data
            else:
                input_buffer = self.fin.read()

        log("Input buffer:", input_buffer)
        if "\n" in input_buffer:
            cmd, rest = input_buffer.split("\n", 1) # split on the first \n only
            self.buffer = rest
        else:
            cmd = input_buffer
            self.buffer = ""

        cmd = cmd.rstrip()
        log("Reader returned cmd:", cmd)
        return cmd

    def close(self):
        log("Reader closed")
        if not self.blocking:
            self.fin.close()

def shutdown():
    ## There may be races since this is called through the signal handling
    log("PaSh daemon is shutting down...")
    if config.bash_mirror is not None:
        ## TODO: Do we need force
        ret = config.bash_mirror.terminate(force=True)
        
        ## The mirror was terminated successfully
        assert(ret)

        config.bash_mirror.wait()
    log("PaSh daemon shut down successfully...")

def main():
    args = init()
    scheduler = Scheduler()
    scheduler.run(in_filename=args.input, out_filename=args.output)
   

if __name__ == "__main__":
    main()
