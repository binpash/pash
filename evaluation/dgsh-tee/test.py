import sys
from subprocess import PIPE, run
import os
GIT_TOP_CMD = [ 'git', 'rev-parse', '--show-toplevel', '--show-superproject-working-tree']
if 'PASH_TOP' in os.environ:
    PASH_TOP = os.environ['PASH_TOP']
else:
    PASH_TOP = run(GIT_TOP_CMD, stdout=PIPE, stderr=PIPE, universal_newlines=True).stdout.rstrip()
sys.path.append(PASH_TOP)

TESTFILES = [f"{PASH_TOP}/evaluation/scripts/input/1M.txt", f"{PASH_TOP}/evaluation/scripts/input/10M.txt", f"{PASH_TOP}/evaluation/scripts/input/100M.txt", f"{PASH_TOP}/evaluation/scripts/input/1G.txt", f"{PASH_TOP}/evaluation/scripts/input/3G.txt", f"{PASH_TOP}/evaluation/scripts/input/100G.txt"]
BATCHSZ = [10000, 100000, 1000000, 10000000, 30000000, 100000000]


from scripts.test_eval.tester import Tests

'''
Make sure to compile dgsh-tee with -DDEBUG 
The test hangs for some reason otherwise.
'''
import os
def run_tests():
    test10M = Tests(TESTFILES[3], BATCHSZ[3])
    tests = ['minimal_sort', 'topn', 'wf', 'diff', 'set-diff', 'double_sort', 'sort', 'shortest_scripts'] #'bigrams', 'spell'
    # tests = ["double_sort", "diff"]
    test_files = []
    for test in tests:
        test_path = f"{PASH_TOP}/evaluation/microbenchmarks/{test}.sh"
        test_files.append(test_path)
    tmp_folder = "tmp_width4_1G"
    test10M.run_test_list(test_files, r_split=True, width=4, log_folder=tmp_folder, dgsh_tee=True)
    test10M.run_test_list(test_files, width=4, log_folder=tmp_folder, dgsh_tee=True)
    test10M.run_test_list(test_files, r_split=True, width=4, log_folder=tmp_folder, dgsh_tee=False)
    test10M.run_test_list(test_files, width=4, log_folder=tmp_folder, dgsh_tee=False)
    print(test10M.get_df()[["test_name", "no_eager", "split_type", "exec_time", "cpu%", "width", "dgsh_tee"]].to_string(index = False))

if __name__ == '__main__':
    os.environ["PASH_TOP"] = PASH_TOP
    run_tests()