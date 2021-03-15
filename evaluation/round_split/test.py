from subprocess import PIPE, run
from resource import getrusage as resource_usage, RUSAGE_SELF
import argparse
import os
from time import time as timestamp
from functools import wraps
import pandas as pd

TESTFILES = ["../scripts/input/1M.txt", "../scripts/input/10M.txt", "../scripts/input/100M.txt", "../scripts/input/1G.txt"]
BATCHSZ = [10000, 100000, 1000000, 10000000]

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
    pash_command = ["../../pa.sh", filename, "-d 1"]
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
    
    def get_dataframe(self, testname, r_data, sys_data, pash_data):
        r_data["test"] = "r-split " + testname
        sys_data["test"] = "system " + testname
        pash_data["test"] = "pash " + testname

        self.df = self.df.append(r_data, ignore_index=True)
        self.df = self.df.append(sys_data, ignore_index=True)
        self.df = self.df.append(pash_data, ignore_index=True)

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
        #round-split
        command = ["./r-wc.sh", self.test_file, self.batch_sz]
        out, data = self.time(command)

        os.environ["IN"] = self.test_file
        #default
        command = ["./wc.sh"]
        out1, data1 = self.time(command)
        # print(out)
        # print(out.split(), out1.split())
        # add assert when merge works

        #pa.sh
        os.environ["IN"] = self.test_file
        pash_file = get_pash_script("wc.sh")
        command = [pash_file]
        out2, data2 = self.time(command)
 
        # assert(out == out1 == out2)
        
        df = self.get_dataframe("wc", data, data1, data2)
        
        return df

    def minimal_grep(self):
        #round-split
        command = ["./r-minimal_grep.sh", self.test_file, self.batch_sz]
        out, data = self.time(command)

        os.environ["IN"] = self.test_file
        #default
        command = ["./minimal_grep.sh"]
        out1, data1 = self.time(command)
        assert(out == out1)

        # pa.sh
        pash_file = get_pash_script("./minimal_grep.sh")
        command = [pash_file]
        out2, data2 = self.time(command)
        df = self.get_dataframe("minimal_grep", data, data1, data2)
        
        return df

    def bell_grep(self):
        #round-split
        command = ["./r-bell_grep.sh", self.test_file, self.batch_sz]
        out, data = self.time(command)

        os.environ["IN"] = self.test_file
        #default
        command = ["./bell_grep.sh"]
        out1, data1 = self.time(command)
        
        assert(out == out1)
        # pa.sh
        pash_file = get_pash_script("bell_grep.sh")
        command = [pash_file]
        out2, data2 = self.time(command)
        df = self.get_dataframe("bell", data, data1, data2)
        
        return df

    def sort(self):
        #round-split
        command = ["./r-sort.sh", self.test_file, self.batch_sz]
        out, data = self.time(command)

        os.environ["IN"] = self.test_file
        #default
        command = ["../microbenchmarks/sort.sh"]
        out1, data1 = self.time(command)

        assert(out == out1)
        # pa.sh
        pash_file = get_pash_script("../microbenchmarks/sort.sh")
        command = [pash_file]
        out2, data2 = self.time(command)
        df = self.get_dataframe("sort", data, data1, data2)
        
        return df

def run_tests():
    print("-----------Running 100M tests---------------")
    test100M = Tests(TESTFILES[2], BATCHSZ[2])
    test100M.bell_grep()
    test100M.sort()
    test100M.wcTest()
    print(test100M.df[["test", "real", "user", "sys", "cpu%"]].to_string(index = False))

    print("\n-----------Running 1G tests---------------")
    test1G = Tests(TESTFILES[3], BATCHSZ[3])
    test1G.bell_grep()
    test1G.sort()
    test1G.wcTest()
    print(test1G.df[["test", "real", "user", "sys", "cpu%"]].to_string(index = False))

    print("\n-----------Running 10M tests---------------")
    test10M = Tests(TESTFILES[1], BATCHSZ[1])
    test10M.bell_grep()
    test10M.sort()
    test10M.wcTest()
    test10M.minimal_grep()
    print(test10M.df[["test", "real", "user", "sys", "cpu%"]].to_string(index = False))

if __name__ == '__main__':
    run_tests()