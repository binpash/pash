from .logparser import LogParser, DEFAULT_LOG_FOLDER
from subprocess import PIPE, run
import argparse
import os
import pandas as pd
import uuid

GIT_TOP_CMD = [
    "git",
    "rev-parse",
    "--show-toplevel",
    "--show-superproject-working-tree",
]
if "PASH_TOP" in os.environ:
    PASH_TOP = os.environ["PASH_TOP"]
else:
    PASH_TOP = run(
        GIT_TOP_CMD, stdout=PIPE, stderr=PIPE, universal_newlines=True
    ).stdout.rstrip()


class Tests(LogParser):
    def __init__(self, in_file=None, batch_sz=100000):
        self.in_file = in_file
        self.batch_sz = str(batch_sz)
        self.log_parser = LogParser()

    def time(self, command, env, stdout=PIPE):
        time_command = ["/usr/bin/time", "-v", "bash"]
        time_command.extend(command)
        result = run(
            time_command,
            stdout=PIPE,
            universal_newlines=True,
            stdin=None,
            stderr=PIPE,
            env=env,
        )
        return result

    def get_df(self):
        return self.log_parser.get_df()

    def run_test(
        self,
        test_path,
        width=2,
        r_split=False,
        batch_size=None,
        in_file=None,
        no_eager=False,
        dgsh_tee=False,
        log_folder=DEFAULT_LOG_FOLDER,
    ):
        if in_file == None:
            in_file = self.in_file

        new_env = os.environ.copy()
        if in_file == None:
            in_file = self.in_file

        new_env["IN"] = in_file
        new_env["PASH_TOP"] = PASH_TOP

        command = [
            f"{PASH_TOP}/pa.sh",
            test_path,
            "--output_time",
            f"-w {width}",
            "-d 1",
        ]

        if r_split:
            command.append("--r_split")
            batch_size = (
                str(batch_size) if batch_size else self.batch_sz
            )  # str(int(os.path.getsize(in_file)/90))
            command.append("--r_split_batch_size")
            command.append(batch_size)
        if no_eager:
            command.append("--no_eager")
        if dgsh_tee:
            command.append("--dgsh_tee")

        result = self.time(command, new_env)

        # add IN file to log
        result.stderr = f"IN {in_file}\n" + result.stderr

        # write stderr to log_file if provided
        log_file = self.__get_log_file__(test_path, log_folder)
        with open(log_file, "w") as f:
            f.write(result.stderr)

        if result.returncode != 0:
            print(f"failed running: {test_path}")
            if log_file:
                print(f"log in {log_file}")

        df = self.log_parser.parse_log(result.stderr)

        return result, df

    # Run provided tests in folder x with the env files
    def run_folder_tests(
        self,
        tests,
        folder,
        width=2,
        r_split=False,
        batch_size=None,
        in_file=None,
        no_eager=False,
        dgsh_tee=False,
        log_folder=None,
    ):
        pass

    # run a list of tests, each test should be the full path of .sh file
    # if log_folder provided it generates unique name for each log
    def run_test_list(
        self,
        tests,
        width=2,
        r_split=False,
        batch_size=None,
        in_file=None,
        no_eager=False,
        dgsh_tee=False,
        log_folder=DEFAULT_LOG_FOLDER,
    ):
        df = pd.DataFrame()
        for test in tests:
            result, dfnew = self.run_test(
                test,
                width,
                r_split,
                batch_size,
                in_file,
                no_eager,
                dgsh_tee,
                log_folder,
            )
            df = df.append(dfnew, ignore_index=True)

        return df

    def __get_log_file__(self, test_path, log_folder):
        if not os.path.exists(log_folder):
            os.makedirs(log_folder, exist_ok=True)

        temp_filename = (
            os.path.basename(test_path).replace(".sh", "")
            + "_"
            + str(uuid.uuid4())
            + ".log"
        )
        log_file = os.path.join(log_folder, temp_filename)
        return log_file
