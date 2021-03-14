from log_parser import LogParser, DEFAULT_LOG_FOLDER
from subprocess import PIPE, run
import argparse
import os
import pandas as pd
import uuid

PASH_TOP = "/home/tamlu/pash"
TESTFILES = [f"{PASH_TOP}/evaluation/scripts/input/1M.txt", f"{PASH_TOP}/evaluation/scripts/input/10M.txt", f"{PASH_TOP}/evaluation/scripts/input/100M.txt", f"{PASH_TOP}/evaluation/scripts/input/1G.txt"]
BATCHSZ = [10000, 100000, 1000000, 10000000]

class Tests(LogParser):
    def __init__(self, in_file = TESTFILES[1], batch_sz = 100000):
        self.in_file = in_file
        self.batch_sz = str(batch_sz)
        self.log_parser = LogParser()

    def time(self, command, env):
        time_command = ["/usr/bin/time" , "-v", "bash"]
        time_command.extend(command)

        result = run(time_command, stdout=PIPE, universal_newlines=True, stdin=None, stderr=PIPE, env=env)
 
        return result

    def get_df(self):
        return self.log_parser.get_df()

    def run_test(self, test_path, width = 2, r_split=False, batch_size=None, in_file=None, no_eager=False, log_folder=DEFAULT_LOG_FOLDER):
        if in_file==None:
            in_file = self.in_file
        
        new_env = os.environ.copy()
        if in_file == None:
            in_file = self.in_file

        new_env["IN"] = in_file

        command = [f"{PASH_TOP}/pa.sh", test_path, "--output_time", f"-w {width}", "-d 1"]

        if r_split:
            command.append("--r_split")
            batch_size = str(batch_size) if batch_size else self.batch_sz #str(int(os.path.getsize(in_file)/90))
            command.append("--r_split_batch_size")
            command.append(batch_size)
        if no_eager:
            command.append("--no_eager")

        result = self.time(command, new_env)
        
        #add IN file to log
        result.stderr = f"IN {in_file}\n" + result.stderr

        #write stderr to log_file if provided
        log_file = self.__get_log_file__(test_path, log_folder)
        with open(log_file, 'w') as f:
            f.write(result.stderr)
            

        if result.returncode != 0:
            print(f"failed running: {test_path}")
            if log_file:
                print(f"log in {log_file}")
        
        df = self.log_parser.parse_log(result.stderr)
        
        return result, df

    #Run provided tests in folder x with the env files
    def run_folder_tests(self, tests, folder, width = 2, r_split=False, batch_size=None, in_file=None, no_eager=False, log_folder=None):
        pass

    #run a list of tests, each test should be the full path of .sh file
    #if log_folder provided it generates unique name for each log
    def run_test_list(self, tests, width = 2, r_split=False, batch_size=None, in_file=None, no_eager=False, log_folder=DEFAULT_LOG_FOLDER):
        df = pd.DataFrame()
        for test in tests:
            result, dfnew = self.run_test(test, width, r_split, batch_size, in_file, no_eager, log_folder)
            df = df.append(dfnew, ignore_index=True)

        return df
    
    def __get_log_file__(self, test_path, log_folder):
        if not os.path.exists(log_folder):
            os.makedirs(log_folder, exist_ok=True)

        temp_filename = os.path.basename(test_path).replace(".sh", "") + "_" + str(uuid.uuid4()) + ".log"
        log_file = os.path.join(log_folder, temp_filename)
        return log_file


def run_tests():
    test10M = Tests(TESTFILES[3], BATCHSZ[3])
    tests = ['minimal_sort', 'topn', 'wf', 'diff', 'set-diff', 'double_sort', 'sort'] #'bigrams', 'spell' 'shortest_scripts'
    
    test_files = []
    for test in tests:
        test_path = f"{PASH_TOP}/evaluation/microbenchmarks/{test}.sh"
        test_files.append(test_path)
    tmp_folder = "tmp_1G_width4/"
    test10M.run_test_list(test_files, r_split=True, width=4, log_folder=tmp_folder)
    test10M.run_test_list(test_files, r_split=True, width=4, no_eager=True, log_folder=tmp_folder)
    test10M.run_test_list(test_files, width=4, log_folder=tmp_folder)
    test10M.run_test_list(test_files, width=4, no_eager=True, log_folder=tmp_folder)
    print(test10M.get_df()[["test_name", "no_eager", "split_type", "exec_time", "cpu%", "width"]].to_string(index = False))

if __name__ == '__main__':
    os.environ["PASH_TOP"] = PASH_TOP
    run_tests()