from log_parser import LogParser
from subprocess import PIPE, run
import argparse
import os
import pandas as pd

PASH_TOP = "/home/tamlu/pash"
TESTFILES = [f"{PASH_TOP}/evaluation/scripts/input/1M.txt", f"{PASH_TOP}/evaluation/scripts/input/10M.txt", f"{PASH_TOP}/evaluation/scripts/input/100M.txt", f"{PASH_TOP}/evaluation/scripts/input/1G.txt"]
BATCHSZ = [10000, 100000, 1000000, 10000000]

class Tests(LogParser):
    def __init__(self, in_file = TESTFILES[1], batch_sz = 100000, df = None):
        self.in_file = in_file
        self.batch_sz = str(batch_sz)
        self.log_parser = LogParser()

    def time(self, command, env):
        time_command = ["/usr/bin/time" , "-v", "bash"]
        time_command.extend(command)

        result = run(time_command, stdout=PIPE, universal_newlines=True, stdin=None, stderr=PIPE, env=env)

        if result.returncode != 0:
            print(result.stderr)
        
        df = self.log_parser.parse_log(result.stderr)
 
        return result.stdout, df

    def get_df(self):
        return self.log_parser.get_df()

    def run_test(self, test_path, width = 2, r_split=False, batch_size=None, in_file=None, no_eager=False, persist_log=None):
        if in_file==None:
            in_file = self.in_file
        
        new_env = os.environ.copy()
        new_env["IN"] = in_file if in_file else self.in_file

        command = [f"{PASH_TOP}/pa.sh", test_path, "--output_time", f"-w {width}", "-d 1"]

        if r_split:
            command.append("--r_split")
            batch_size = batch_size if batch_size else self.batch_sz
            command.append("--r_split_batch_size")
            command.append(batch_size)
        if no_eager:
            command.append("--no_eager")

        out, df = self.time(command, new_env)
        return out, df

    #not recommended for now, Add timeout?
    def run_folder_tests(self, test_folder, width = 2, r_split=False, batch_size=None, in_file=None, no_eager=False, persist_log=None):
        pass

    #run a list of tests, each test should be the full path of .sh file
    def run_test_list(self, tests, width = 2, r_split=False, batch_size=None, in_file=None, no_eager=False, persist_log=None):
        for test in tests:
            out, df = self.run_test(test, width, r_split, batch_size, in_file, no_eager, persist_log)
        
def run_tests():
    test10M = Tests(TESTFILES[1], BATCHSZ[1])
    tests = ['minimal_grep', 'minimal_sort', 'topn', 'wf', 'diff', 'set-diff', 'double_sort', 'shortest_scripts'] #'bigrams', 'spell'
    
    test_files = []
    for test in tests:
        test_path = f"{PASH_TOP}/evaluation/microbenchmarks/{test}.sh"
        test_files.append(test_path)
    
    test10M.run_test_list(test_files, r_split=True, width=4)
    test10M.run_test_list(test_files, width=4)
    
    print(test10M.get_df()[["test_name", "no_eager", "split_type", "exec_time", "cpu%", "width"]].to_string(index = False))

if __name__ == '__main__':
    os.environ["PASH_TOP"] = PASH_TOP
    run_tests()