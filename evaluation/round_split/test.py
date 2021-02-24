from subprocess import PIPE, run
from resource import getrusage as resource_usage, RUSAGE_SELF
import argparse
import os
from time import time as timestamp
from functools import wraps

TESTFILES = ["/home/ubuntu/pash/evaluation/scripts/input/1M.txt", "/home/ubuntu/pash/evaluation/scripts/input/10M.txt", "/home/ubuntu/pash/evaluation/scripts/input/100M.txt", "/home/ubuntu/pash/evaluation/scripts/input/1G.txt"]


def process_gnu_time(time_data):
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


class Tests():
    ARGS = ""
    def __init__(self, test_file = TESTFILES[1], batch_sz = 100000):
        self.test_file = test_file
        self.batch_sz = str(batch_sz)
    
    def time(self, command):
        time_command = ["/usr/bin/time" , "-v"]
        time_command.extend(command)
        start_time = timestamp()
        result = run(time_command, stdout=PIPE, universal_newlines=True, stdin=None, stderr=PIPE)
        end_time = timestamp()
        data = process_gnu_time(result.stderr)
        data["real"] = end_time - start_time
        return result.stdout, data

    def wcTest(self):
        command = ["./r-wc.sh", self.test_file, self.batch_sz]
        out, data = self.time(command)

        command = ["wc", self.test_file]
        out1, data1 = self.time(command)

        return data
    


if __name__ == '__main__':
    tests = Tests()
    tests.wcTest()