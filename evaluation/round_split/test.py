import sys
from subprocess import PIPE, run
import os
GIT_TOP_CMD = [ 'git', 'rev-parse', '--show-toplevel', '--show-superproject-working-tree']
if 'PASH_TOP' in os.environ:
    PASH_TOP = os.environ['PASH_TOP']
else:
    PASH_TOP = run(GIT_TOP_CMD, stdout=PIPE, stderr=PIPE, universal_newlines=True).stdout.rstrip()
sys.path.append(PASH_TOP)

TEST_FILE_DIR = f"{PASH_TOP}/evaluation/benchmarks/oneliners/input"
TESTFILES = [f"{TEST_FILE_DIR}/1M.txt", f"{TEST_FILE_DIR}/10M.txt", f"{TEST_FILE_DIR}/100M.txt", f"{TEST_FILE_DIR}/1G.txt",f"{TEST_FILE_DIR}/3G.txt"]
BATCHSZ = [32000, 100000, 1000000, 10000000]


from scripts.test_eval.tester import Tests
from scripts.test_eval.logparser import LogParser
import os


def get_tests(dir):
    tests = []
    for test in os.listdir(dir):
        if test.endswith(".sh"):
            tests.append(os.path.join(dir, test))
    return tests

# def run_bash(dir, in_file, save_logs_to="tmp_bash/"):
#     tests = get_tests(dir)
#     tester = Tests(in_file)
#     # command = [tmp]
    


def run_tests(dir):
    tester = Tests(TESTFILES[3], BATCHSZ[2])
    test_files = get_tests(dir)
    # tests = ['minimal_sort', 'topn', 'wf', 'diff', 'set-diff', 'double_sort', 'sort'] #'bigrams', 'spell' 'shortest_scripts'
    # # tests = ["double_sort", "diff"]
    for width in [4, 8]:
        batch_size = BATCHSZ[2]
        tmp_folder = f"tmp_width{width}_1G_{batch_size}bz_nospill_5Mbuffer/"
        df = tester.run_test_list(test_files, r_split=True, width=width, log_folder=tmp_folder, batch_size=batch_size, dgsh_tee=True)
        print(f"-------Width {width} BatchSize {batch_size}")
        print(df[["test_name", "no_eager", "split_type", "exec_time", "cpu%", "width"]].to_string(index = False))
        print("\n\n")
    print(tester.get_df()[["test_name", "no_eager", "split_type", "exec_time", "cpu%", "width"]].to_string(index = False))

if __name__ == '__main__':
    os.environ["PASH_TOP"] = PASH_TOP
    run_tests(f"{PASH_TOP}/evaluation/round_split/scripts/")
    # run_bash(f"{PASH_TOP}/evaluation/round_split/scripts/")