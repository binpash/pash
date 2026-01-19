import pandas as pd
import argparse
import re
import os
from sys import argv

DEFAULT_LOG_FOLDER = "tmp_log/"


class LogParser:
    """
    A class used to parse the pa.sh log files

    All parse_* methods return a dataframe of only the files parsed in this call.
    Use get_df for all parsed files across multible calls to parse_*.

    Methods:
        parse_file: parses a log file
        parse_folder: parses log files in a folder
        parse_log: parses a given log string
        get_df: returns a comprehensive dataframe of every
                log parsed (using any of the functions above)
                during the function lifetime.

    Dataframe columns:
        - test_name
        - IN
        - split_type
        - no_eager
        - dgsh_tee
        - width
        - r_split_batch_size
        - exec_time
        - backend_time
        - execiton_time
        - compilation_time
        - preprocess_time
        - eager_nodes
        - compiler_exit
        - gnu_real
        - gnu_usr
        - gnu_sys
        - cpu%
        - max_res_size
        - avg_res_size
        - major_pagefaults
        - minor_pagefaults
    """

    def __init__(self, df=None):
        self.df = df if df else pd.DataFrame()

    def parse_log(self, log: str) -> pd.DataFrame:
        """
        Parses a pa.sh log with path file_path
        Return:
            A single entry pandas dataframe, or None if failed
        """

        border = "-" * 40
        argslog, pashlog, timelog = log.split(border)

        args_of_interest = set(
            [
                "input",
                "width",
                "output_time",
                "no_eager",
                "r_split",
                "r_split_batch_size",
                "IN",
                "dgsh_tee",
            ]
        )
        parsed_args = self.__parse_args__(argslog, args_of_interest)

        tags_of_interest = set(
            [
                "Execution time",
                "Backend time",
                "Compilation time",
                "Preprocessing time",
                "Eager nodes",
                "Compiler exited with code",
            ]
        )
        parsed_log = self.__parse_pash_log__(pashlog, tags_of_interest)

        # can be empty
        parsed_time = self.__parse_time_log__(timelog)

        if not parsed_args["input"]:
            return None

        df = pd.DataFrame()

        test_name = os.path.basename(parsed_args["input"]).replace(".sh", "")
        split_type = "r-split" if parsed_args["r_split"] else "auto-split"

        data = {
            # From Args
            "test_name": test_name,
            "IN": os.path.basename(parsed_args["IN"]),
            "split_type": split_type,
            "no_eager": parsed_args["no_eager"],
            "width": int(parsed_args["width"]),
            "r_split_batch_size": int(parsed_args["r_split_batch_size"]),
            "dgsh_tee": parsed_args["dgsh_tee"],
            # From pash log
            "exec_time": parsed_log["Execution time"],
            "backend_time": parsed_log["Backend time"],
            "compilation_time": parsed_log["Compilation time"],
            "preprocess_time": parsed_log["Preprocessing time"],
            "eager_nodes": int(parsed_log["Eager nodes"]),
            "compiler_exit": parsed_log["Compiler exited with code"],
            # From time
            "gnu_real": parsed_time["gnu_real"],
            "gnu_usr": parsed_time["user"],
            "gnu_sys": parsed_time["sys"],
            "cpu%": parsed_time["cpu%"],
            "max_res_size": int(parsed_time["max_resident"]),
            "avg_res_size": int(parsed_time["average_resident"]),
            "major_pagefaults": int(parsed_time["major_pagefaults"]),
            "minor_pagefaults": int(parsed_time["minor_pagefaults"]),
        }

        # update local and global df
        df = df.append(data, ignore_index=True)
        self.df = self.df.append(data, ignore_index=True)

        return df

    def parse_file(self, log_file: str) -> pd.DataFrame:
        """
        Parses a pa.sh log with path file_path
        Return:
            A single entry pandas dataframe
        """
        try:
            with open(log_file, "r") as f:
                log = f.read()
                df = self.parse_log(log)
                return df
        except:
            print("failed to parse", log_file)
            return pd.DataFrame()

    def parse_folder(self, path: str) -> pd.DataFrame:
        """
        Parses all valid files ending with .log in the path directory.
        Params:
            path: should be a path to a directory containing .log files
        Return:
            pandas dataframe with all parsed logs
        """
        log_files = [
            os.path.join(path, f) for f in os.listdir(path) if f.endswith(".log")
        ]
        ret_df = pd.DataFrame()
        for log_file in log_files:
            df = self.parse_file(log_file)
            ret_df = ret_df.append(df, ignore_index=True)

        return ret_df

    def get_df(self):
        self.df["no_eager"] = self.df["no_eager"].astype(bool)
        self.df["width"] = self.df["width"].astype(int)
        self.df["r_split_batch_size"] = self.df["r_split_batch_size"].astype(int)
        self.df["dgsh_tee"] = self.df["dgsh_tee"].astype(bool)
        return self.df

    def __parse_args__(self, args: str, args_of_interest):
        lines = args.split("\n")
        args_dict = {i: False for i in args_of_interest}
        for line in lines:
            try:
                arg, val = line.split(" ")
                if arg in args_of_interest:
                    if val == "True" or val == "False":
                        args_dict[arg] = True if val == "True" else False
                    else:
                        args_dict[arg] = val
            except:
                continue
        return args_dict

    def __parse_pash_log__(self, args: str, tags_of_interest):
        lines = args.split("\n")
        log_dict = {i: 0 for i in tags_of_interest}

        for line in lines:
            try:
                tag, val = line.split(": ")
                if tag in tags_of_interest:
                    if tag.endswith("time"):
                        log_dict[tag] += float(val.split(" ")[0])
                    else:
                        log_dict[tag] = val
            except:
                continue
        return log_dict

    def __parse_time_log__(self, timelog: str):
        data_start = timelog.find("Command being timed: ")
        time_data = timelog[data_start:]

        lines = [
            line.split(": ")[1] for line in time_data.replace("\t", "").split("\n")[:-1]
        ]
        if len(lines) < 23:
            lines = [False] * 23
        data = {
            "command": lines[0],
            "user": lines[1],
            "sys": lines[2],
            "cpu%": lines[3],
            "gnu_real": lines[4],
            "max_resident": lines[9],
            "average_resident": lines[10],
            "major_pagefaults": lines[11],
            "minor_pagefaults": lines[12],
            "exit_status": lines[22],
        }
        return data


# can be used in case we only can parse the time (default commands)
def process_gnu_time(time_data):
    data_start = time_data.find("Command being timed: ")
    time_data = time_data[data_start:]
    lines = [
        line.split(": ")[1] for line in time_data.replace("\t", "").split("\n")[:-1]
    ]
    data = {
        "command": lines[0],
        "user": lines[1],
        "sys": lines[2],
        "cpu%": lines[3],
        "gnu_real": lines[4],
        "max_resident": lines[9],
        "average_resident": lines[10],
        "major_pagefault": lines[11],
        "minor_pagefault": lines[12],
        "exit_status": lines[22],
    }
    return data


if __name__ == "__main__":
    # sample execution
    log_parser = LogParser()
    # can pass folder name in first argument
    if len(argv) > 1:
        df = log_parser.parse_folder(argv[1])
    else:
        df = log_parser.parse_folder(DEFAULT_LOG_FOLDER)
    print(
        log_parser.get_df()[
            [
                "test_name",
                "IN",
                "r_split_batch_size",
                "no_eager",
                "split_type",
                "exec_time",
                "cpu%",
                "width",
            ]
        ].to_string(index=False)
    )
