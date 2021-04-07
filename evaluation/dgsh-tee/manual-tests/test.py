from subprocess import PIPE, run
from resource import getrusage as resource_usage, RUSAGE_SELF
import argparse
import os
from time import time as timestamp
from functools import wraps
import pandas as pd

PASH_TOP = os.environ['PASH_TOP']
TESTFILES = [f"{PASH_TOP}/evaluation/scripts/input/1M.txt", f"{PASH_TOP}/evaluation/scripts/input/10M.txt", f"{PASH_TOP}/evaluation/scripts/input/100M.txt", f"{PASH_TOP}/evaluation/scripts/input/1G.txt", f"{PASH_TOP}/evaluation/scripts/input/3G.txt", f"{PASH_TOP}/evaluation/scripts/input/100G.txt"]
BATCHSZ = [10000, 100000, 1000000, 10000000, 30000000, 100000000]

def process_gnu_time(time_data):
    data_start = time_data.find("Command being timed: ")
    time_data = time_data[data_start: ]
    lines = [line.split(": ")[1] for line in time_data.replace("\t", "").split("\n")[:-1]]
    data = {
        "command" : lines[0],
        "user" : lines[1],
        "sys" : lines[2],
        "cpu%" : lines[3],
        "gnu_real" : lines[4],
        "max_resident" : lines[9],
        "average_resident": lines[10],
        "major_pagefault" : lines[11],
        "minor_pagefault" : lines[12],
        "exit_status" : lines[22]
    }
    return data

def get_pash_script(filename) :
    pash_command = [f"{PASH_TOP}/pa.sh", filename, "-d 1"]
    result = run(pash_command, stdout=PIPE, universal_newlines=True, stdin=None, stderr=PIPE)
    data = result.stderr.split("\n")
    for line in data:
        if "Optimized script saved in:" in line:
            script = line.split(": ")[1]
            # run(["chmod", "+x", script])
            return script
    print("error getting pash optimized script")
    return None


class Tests():
    ARGS = ""
    def __init__(self, test_file = TESTFILES[1], batch_sz = 100000, df = None):
        self.test_file = test_file
        self.batch_sz = str(batch_sz)
        self.df = df if df else pd.DataFrame()
    
    def get_dataframe(self, testname, all_data):
        for test in all_data:
            test_data = all_data[test]
            test_data["test"] = test + " " + testname
            self.df = self.df.append(test_data, ignore_index=True)
        return self.df

    def time(self, command):
        time_command = ["/usr/bin/time" , "-v", "bash"]
        time_command.extend(command)

        start_time = timestamp()
        result = run(time_command, stdout=PIPE, universal_newlines=True, stdin=None, stderr=PIPE)
        end_time = timestamp()
        if result.returncode != 0:
            print(result.stderr)
        data = process_gnu_time(result.stderr)
        data["real"] = end_time - start_time
        return result.stdout, data

    
    def wcTest(self):
        data = {}
        #round-split
        command = ["./dgsh-wc.sh", self.test_file, self.batch_sz]
        out, data["dgsh-split"] = self.time(command)

        os.environ["IN"] = self.test_file
        #default
        command = [f"{PASH_TOP}/evaluation/round_split/manual-tests/wc.sh"]
        out1, data["system"] = self.time(command)

        # TODO: add assert when merge works

        #pa.sh
        os.environ["IN"] = self.test_file
        pash_file = get_pash_script(f"{PASH_TOP}/evaluation/round_split/manual-tests/wc.sh")
        command = [pash_file]
        out2, data["pash"] = self.time(command)
 
        # assert(out == out1 == out2)
        
        df = self.get_dataframe("wc", data)
        
        return df


    def sort(self):
        data = {}
        #round-split
        command = ["./dgsh-sort.sh", self.test_file, self.batch_sz]
        out, data["dgsh-split"] = self.time(command)

        #using -r
        command = ["./dgsh-raw-sort.sh", self.test_file, self.batch_sz]
        out1, data["dgsh-raw-split"] = self.time(command)

        os.environ["IN"] = self.test_file
        #default
        command = [f"{PASH_TOP}/evaluation/microbenchmarks/sort.sh"]
        out2, data["system"] = self.time(command)

        assert(out == out2)
        assert(out1 == out2)

        # pa.sh
        pash_file = get_pash_script(f"{PASH_TOP}/evaluation/microbenchmarks/sort.sh")
        command = [pash_file]
        out2, data["pash"] = self.time(command)
        df = self.get_dataframe("sort", data)
        
        return df
        
        
def run_tests():
    # print("-----------Running 100M tests---------------")
    # test100M = Tests(TESTFILES[2], BATCHSZ[2])
    # test100M.sort()
    # test100M.wcTest()
    # print(test100M.df[["test", "real", "user", "sys", "cpu%"]].to_string(index = False))

    # print("\n-----------Running 1G tests---------------")
    # test1G = Tests(TESTFILES[3], BATCHSZ[3])
    # test1G.sort()
    # test1G.wcTest()
    # print(test1G.df[["test", "real", "user", "sys", "cpu%"]].to_string(index = False))

    print("\n-----------Running 10M tests---------------")
    test10M = Tests(TESTFILES[1], BATCHSZ[1])
    test10M.sort()
    test10M.wcTest()
    print(test10M.df[["test", "real", "user", "sys", "cpu%"]].to_string(index = False))

if __name__ == '__main__':
    os.environ["PASH_TOP"] = PASH_TOP
    run_tests()