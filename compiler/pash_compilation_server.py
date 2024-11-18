import signal
from threading import Thread
from datetime import datetime, timedelta

# import queue

from sh_expand import env_vars_util

import config
from pash_graphviz import maybe_generate_graphviz
import pash_compiler
from util import *
from dspash.worker_manager import WorkersManager
import server_util

from cli import BaseParser
from custom_error import * 

##
## A Daemon (not with the strict Unix sense)
## that responds to requests for compilation
##


def handler(signum, frame):
    log("Signal:", signum, "caught")
    shutdown()


signal.signal(signal.SIGTERM, handler)


def parse_args():
    parser = BaseParser(add_help=False)
    parser.add_pash_args()
    args, _ = parser.parse_known_args()

    return args


# Initialize the daemon


def init():
    ## Set the logging prefix
    config.LOGGING_PREFIX = "Daemon: "

    args = parse_args()
    config.set_config_globals_from_pash_args(args)

    # Load the configuration
    if not config.config:
        config.load_config(args.config_path)

    pash_compiler.runtime_config = config.config["distr_planner"]

    return args


##
## This class holds information for each process id
##
class ProcIdInfo:
    def __init__(self, input_ir, compiler_config, exec_time=None, start_exec_time=None):
        self.input_ir = input_ir
        self.compiler_config = compiler_config
        self.exec_time = exec_time
        self.start_exec_time = start_exec_time
        ## TODO: Extend it with other info from scheduler, like dependencies

    def set_exec_time(self, exec_time):
        self.exec_time = exec_time

    def set_start_exec_time(self, start_exec_time):
        self.start_exec_time = start_exec_time

    def get_start_exec_time(self):
        return self.start_exec_time

    def __repr__(self):
        return f"ProcIdInfo(InputIR:{self.input_ir}, CompConfig:{self.compiler_config}, ExecTime:{self.exec_time})"


class Scheduler:
    """Takes care of running processes in parallel if there is no conflict.
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
        self.process_resources = (
            {}
        )  # map process_id -> (input_resources, output_resources)
        self.next_id = 0
        self.running_procs = 0
        self.unsafe_running = False
        self.done = False
        self.cmd_buffer = ""
        self.connection_manager = None
        self.reader_pipes_are_blocking = True
        self.request_processing_start_time = 0
        ## TODO: Make that be a class or something

        ## A map that keeps mappings between proc_id and (input_ir, width, exec_time)
        self.process_id_input_ir_map = {}
        ## This is a map from input IRs, i.e., locations in the code, to a list of process_ids
        self.input_ir_to_process_id_map = {}

    def check_resources_safety(self, process_id):
        proc_input_resources, proc_output_resources = self.process_resources[process_id]
        all_proc_resources = proc_input_resources.union(proc_output_resources)
        if self.output_resources.intersection(
            all_proc_resources
        ) or self.input_resources.intersection(proc_output_resources):
            return False
        return True

    ##############################################################################
    ##
    ## Profiling-based Compiler-config
    ##
    ##############################################################################

    ## TODO: It is possible that there will be many process ids being gathered here.
    ##       If that becomes an overhead, we need to consider removing old ones
    ##         or keeping sketches of the information (instead of their raw numbers).

    def determine_compiler_config(self, input_ir_file):
        if config.pash_args.profile_driven:
            ## A default value
            selected_width = config.pash_args.width
            # log("Using profiles to select width!")

            ## Goal: Find the highest width that gives benefits
            ##
            ## Strategy, start trying lower widths, if the time seems to drop, keep trying lower.
            ##
            width_avgs = self.get_averages_per_width(input_ir_file)
            log("Width averages:", width_avgs)
            widths = width_avgs.keys()

            ## If we have at least 1, with a specific width,
            ##   and the minimum width has the lowest average, then try one lower
            if len(widths) > 0:
                min_width = min(widths)
                min_width_average = width_avgs[min_width]
                # log("Minimum width:", min_width, "with average:", min_width_average)

                best_width = min_width
                best_width_average = min_width_average
                for width, average in width_avgs.items():
                    ## We find the width with the best average
                    ## If that ends up being the minimum, then we also try one lower.
                    if average < best_width_average:
                        best_width = width
                        best_width_average = average

                if best_width == min_width and min_width > 1:
                    ## Divide the min_width by 2 and try again
                    selected_width = min_width // 2
                    log(
                        "Best width is the lowest width, trying with width:",
                        selected_width,
                    )
                else:
                    selected_width = best_width
                    log("Best width is:", best_width, "We will keep executing with it.")
        else:
            selected_width = config.pash_args.width

        log("Selected width:", selected_width)
        return pash_compiler.CompilerConfig(selected_width)

    def get_averages_per_width(self, input_ir_file):
        ## If we haven't gathered any statistic yet
        if not input_ir_file in self.input_ir_to_process_id_map:
            return {}

        all_proc_ids = self.input_ir_to_process_id_map[input_ir_file]

        ## TODO: If that becomes an overhead, we need to add a new
        ##       data structure to keep the latest averages (sum, len)
        ##       to avoid recomputing every time.
        width_times = {}
        for proc_id in all_proc_ids:
            proc_info = self.process_id_input_ir_map[proc_id]
            width = proc_info.compiler_config.width
            exec_time = proc_info.exec_time

            if exec_time is not None:
                try:
                    width_times[width].append(exec_time)
                except:
                    width_times[width] = [exec_time]

        ## We have gathered all times for each width
        width_avgs = {}
        for width, exec_times in width_times.items():
            width_avgs[width] = sum(exec_times) / len(exec_times)

        return width_avgs

    ## This adds the time measurement, or just removes the entry if there is no exec_time (for space reclamation)
    def handle_time_measurement(self, process_id, exec_time):
        ## 2023-12-08 KK: When in parallel pipelines we receive two exits (when I tried to make it one something got stuck...)
        ##                so this assert is not true
        # assert self.process_id_input_ir_map[process_id].exec_time is None

        ## If we don't have the exec time we do Nothing
        ##
        ## TODO: Consider removing past entries that have no execution time.
        if exec_time is None:
            pass
        else:
            self.process_id_input_ir_map[process_id].set_exec_time(exec_time)

        # log("All measurements:", self.process_id_input_ir_map)

    def add_proc_id_map(self, process_id, input_ir_file, compiler_config):
        assert not process_id in self.process_id_input_ir_map
        self.process_id_input_ir_map[process_id] = ProcIdInfo(
            input_ir_file, compiler_config
        )

        ## Add the mapping from ir to process_id
        try:
            self.input_ir_to_process_id_map[input_ir_file].append(process_id)
        except:
            self.input_ir_to_process_id_map[input_ir_file] = [process_id]
        # log("All mappings from input_ir to proc_id:", self.input_ir_to_process_id_map)

    ##############################################################################

    def compile_and_add(self, compiled_script_file, var_file, input_ir_file):
        process_id = self.get_next_id()
        run_parallel = False
        compile_success = False
        current_region_parallelizable = True 

        variable_reading_start_time = datetime.now()
        # Read any shell variables files if present
        vars_dict = env_vars_util.read_vars_file(var_file, config.BASH_VERSION)
        config.set_vars_file(var_file, vars_dict)

        variable_reading_end_time = datetime.now()
        print_time_delta(
            "Variable Loading", variable_reading_start_time, variable_reading_end_time
        )

        daemon_compile_start_time = datetime.now()
        ## TODO: Make the compiler config based on profiling data
        compiler_config = self.determine_compiler_config(input_ir_file)
        ## Add the process_id -> input_ir mapping
        self.add_proc_id_map(process_id, input_ir_file, compiler_config)

        # check if any general exceptions are caught to report to --assert_compiler_success flag 
        try: 
            ast_or_ir = pash_compiler.compile_ir(
                input_ir_file, compiled_script_file, config.pash_args, compiler_config
            )
        except NotAllRegionParallelizableError: 
            ast_or_ir = None 
            current_region_parallelizable = False
                    

        daemon_compile_end_time = datetime.now()
        print_time_delta(
            "Daemon Compile", daemon_compile_start_time, daemon_compile_end_time
        )

        self.wait_unsafe()
        if ast_or_ir != None:
            compile_success = True

            maybe_generate_graphviz(
                ast_or_ir, config.pash_args, name=f"dfg-{process_id}"
            )

            proc_input_resources = set(
                map(
                    lambda out: str(out.resource)
                    if str(out.resource) != "None"
                    else out,
                    ast_or_ir.all_input_fids(),
                )
            )
            proc_output_resources = set(
                map(
                    lambda out: str(out.resource)
                    if str(out.resource) != "None"
                    else out,
                    ast_or_ir.all_output_fids(),
                )
            )

            self.process_resources[process_id] = (
                proc_input_resources,
                proc_output_resources,
            )

            run_parallel = self.check_resources_safety(process_id)
            if run_parallel:
                self.input_resources = self.input_resources.union(proc_input_resources)
                self.output_resources = self.output_resources.union(
                    proc_output_resources
                )

        if not run_parallel:
            ## If we are not running in parallel everything has to finish first before scheduling for execution
            self.wait_for_all()
        else:
            ## Wait if we have more pipelines running than our current limit
            self.wait_until_limit(config.pash_args.parallel_pipelines_limit)
        
        if compile_success:
            response = server_util.success_response(
                f"{process_id} {compiled_script_file} {var_file} {input_ir_file}"
            )
        elif not current_region_parallelizable: 
            # send specified message to say current region is not parallelizable instead of general exception caught
            response = server_util.error_response(f"{process_id} current region is not parallelizable; failed to compile")
            self.unsafe_running = True
        else:
            response = server_util.error_response(f"{process_id} failed to compile")
            self.unsafe_running = True
           

        ## Do not increase the running procs if assert_all_regions_parallelizable is enabled
        ##  and compilation failed, since nothing will run then.
        ## Do not increase when compile is not successful but regions are parallelizable (in the case that general exceptions are caught), 
        ##  nothing will run in this case also 
        if (not compile_success and config.pash_args.assert_all_regions_parallelizable): 
            pass 
        elif (not compile_success and current_region_parallelizable and config.pash_args.assert_compiler_success): 
            pass 
        else:
            self.running_procs += 1

        ## Get the time before we start executing (roughly) to determine how much time this command execution will take
        command_exec_start_time = datetime.now()
        self.process_id_input_ir_map[process_id].set_start_exec_time(
            command_exec_start_time
        )
        return response

    def remove_process(self, process_id):
        log("The following process exited:", process_id)
        if process_id in self.process_resources:
            del self.process_resources[process_id]
            # TODO: Should be improved to not rebuild inputs and outputs from scratch maybe use counters
            self.input_resources = set().union(
                *[
                    input_resources
                    for input_resources, _ in self.process_resources.values()
                ]
            )
            self.output_resources = set().union(
                *[
                    output_resources
                    for _, output_resources in self.process_resources.values()
                ]
            )

        self.running_procs -= 1
        if self.running_procs == 0:
            self.unsafe_running = False

    def get_next_id(self):
        self.next_id += 1
        return self.next_id

    def wait_for_all(self):
        log(
            "Waiting for all processes to finish."
        )
        self.wait_until_limit(1)
        self.unsafe_running = False

    def wait_until_limit(self, limit: int):
        log(
            f"Waiting for less than {limit} processes to be running. There are",
            self.running_procs,
            "processes remaining.",
        )
        while self.running_procs >= limit:
            input_cmd = self.get_input()
            # must be exit command or something is wrong
            if input_cmd.startswith("Exit:"):
                self.handle_exit(input_cmd)
            else:
                raise Exception(f"Command should be exit but it was {input_cmd}")

    def handle_exit(self, input_cmd):
        assert input_cmd.startswith("Exit:")
        process_id = int(input_cmd.split(":")[1])

        ## Get the execution time
        command_finish_exec_time = datetime.now()
        command_start_exec_time = self.process_id_input_ir_map[
            process_id
        ].get_start_exec_time()
        exec_time = (command_finish_exec_time - command_start_exec_time) / timedelta(
            milliseconds=1
        )
        log("Process:", process_id, "exited. Exec time was:", exec_time)
        self.handle_time_measurement(process_id, exec_time)
        self.remove_process(process_id)
        ## Necessary so that Exit doesn't block
        self.close_last_connection()

    def wait_unsafe(self):
        log("Unsafe running:", self.unsafe_running)
        if self.unsafe_running:
            assert self.running_procs == 1
            self.wait_for_all()
            self.unsafe_running = False

    def parse_and_run_cmd(self, input_cmd):
        if input_cmd.startswith("Compile"):
            (
                compiled_script_file,
                var_file,
                input_ir_file,
            ) = self.__parse_compile_command(input_cmd)
            response = self.compile_and_add(
                compiled_script_file, var_file, input_ir_file
            )
            request_processing_end_time = datetime.now()
            print_time_delta(
                "Request handling",
                self.request_processing_start_time,
                request_processing_end_time,
            )
            ## Send output to the specific command
            self.respond(response)
        elif input_cmd.startswith("Exit:"):
            self.handle_exit(input_cmd)
        elif input_cmd.startswith("Done"):
            self.wait_for_all()
            ## We send output to the top level pash process
            ## to signify that we are done.
            self.respond("All finished")
            self.done = True
        elif input_cmd.startswith("Daemon Start") or input_cmd == "":
            ## This happens when pa.sh first connects to daemon to see if it is on
            self.close_last_connection()
        else:
            log(server_util.error_response(f"Error: Unsupported command: {input_cmd}"))
            raise Exception(f"Error: Unsupported command: {input_cmd}")

    ## This method calls the reader to get an input
    def get_input(self):
        return self.connection_manager.get_next_cmd()

    ## This method responds to the last connection and closes it
    def respond(self, message):
        self.connection_manager.respond(message)

    ## This method closes the last connection that we got input from
    def close_last_connection(self):
        self.connection_manager.close_last_connection()

    def __parse_compile_command(self, input):
        try:
            components = input.rstrip().split("|")
            compiled_script_file = components[0].split(":")[1]
            var_file = components[1].split(":")[1]
            input_ir_file = components[2].split(":")[1]
            return compiled_script_file, var_file, input_ir_file
        except:
            raise Exception(f"Parsing failure for line: {input}")

    def run(self):
        ## By default communicate through sockets, except if the user wants to do it through pipes
        if config.pash_args.daemon_communicates_through_unix_pipes:
            in_filename = os.getenv("RUNTIME_IN_FIFO")
            out_filename = os.getenv("RUNTIME_OUT_FIFO")
            self.connection_manager = server_util.UnixPipeReader(
                in_filename, out_filename, self.reader_pipes_are_blocking
            )
        else:
            self.connection_manager = server_util.SocketManager(
                os.getenv("DAEMON_SOCKET")
            )
        while not self.done:
            # Process a single request
            input_cmd = self.get_input()
            self.request_processing_start_time = datetime.now()

            ## Parse the command (potentially also sending a response)
            self.parse_and_run_cmd(input_cmd)

        self.connection_manager.close()
        shutdown()


def shutdown():
    ## There may be races since this is called through the signal handling
    log("PaSh daemon is shutting down...")
    log("PaSh daemon shut down successfully...")


def main():
    args = init()
    if args.distributed_exec:
        worker_manager = WorkersManager()
        config.worker_manager = worker_manager
        worker_manager_thread = Thread(target=worker_manager.run)
        worker_manager_thread.start()

    scheduler = Scheduler()
    scheduler.run()


if __name__ == "__main__":
    main()
