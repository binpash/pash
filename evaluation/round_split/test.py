import sys
from subprocess import PIPE, run
import os
GIT_TOP_CMD = [ 'git', 'rev-parse', '--show-toplevel', '--show-superproject-working-tree']
if 'PASH_TOP' in os.environ:
    PASH_TOP = os.environ['PASH_TOP']
else:
    PASH_TOP = run(GIT_TOP_CMD, stdout=PIPE, stderr=PIPE, universal_newlines=True).stdout.rstrip()
sys.path.append(PASH_TOP)

from scripts.test_eval.tester import Tests
from scripts.test_eval.logparser import LogParser
import os

TEST_FILE_DIR = f"{PASH_TOP}/evaluation/benchmarks/oneliners/input"
TESTFILES = [f"{TEST_FILE_DIR}/1M.txt", f"{TEST_FILE_DIR}/10M.txt", f"{TEST_FILE_DIR}/100M.txt", f"{TEST_FILE_DIR}/1G.txt",f"{TEST_FILE_DIR}/3G.txt"]
BATCHSZ = [32000, 100000, 1000000, 10000000]




def get_tests(dir):
    tests = []
    for test in os.listdir(dir):
        if test.endswith(".sh"):
            tests.append(os.path.join(dir, test))
    return tests


def run_tests():
    tester = Tests(TESTFILES[3], BATCHSZ[2])

    tests_dir = f"{PASH_TOP}/evaluation/benchmarks/oneliners/"   
    tests = ['spell.sh', 'set-diff.sh', 'nfa-regex.sh', 'top-n.sh', 'sort.sh', 'bi-grams.sh', 'sort-sort.sh', 'diff.sh', 'shortest-scripts.sh', 'wf.sh']
    test_files = list(map(lambda t: os.path.join(tests_dir, t), tests))
    for width in [8]:
        tester.run_test_list(test_files, r_split=True, width=width, log_folder="roundrobin_split/")
        tester.run_test_list(test_files, r_split=False, width=width, log_folder="sequential_split/")
        
    print(tester.get_df()[["test_name", "split_type", "gnu_real", "cpu%", "width"]].to_string(index = False))

if __name__ == '__main__':
    os.environ["PASH_TOP"] = PASH_TOP
    run_tests()
