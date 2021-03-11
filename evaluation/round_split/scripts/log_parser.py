import pandas as pd
import argparse
import re
import os

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
        - split_type
        - no_eager
        - width
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
        
    def parse_log(self, log: str)->pd.DataFrame:
        """
        Parses a pa.sh log with path file_path
        Return:
            A single entry pandas dataframe, or None if failed
        """
        
        border = "-"*40
        argslog, pashlog, timelog = log.split(border)

        args_of_interest = set(["input", "width", "output_time", "no_eager", "r_split", "r_split_batch_size"])
        parsed_args = self.__parse_args__(argslog, args_of_interest)

        tags_of_interest = set(["Execution time", "Backend time", "Compilation time", "Preprocessing time", "Eager nodes", "Compiler exited with code"])
        parsed_log = self.__parse_pash_log__(pashlog, tags_of_interest)

        #can be empty
        parsed_time = self.__parse_time_log__(timelog)

        if not parsed_args["input"]:
            return None

        df = pd.DataFrame()

        test_name = os.path.basename(parsed_args["input"]).replace(".sh", "")
        split_type = "r-split" if parsed_args["r_split"]=="True" else "auto-split"

        data = {
            #From Args
            "test_name" : test_name,
            "split_type" : split_type,
            "no_eager" : parsed_args["no_eager"],
            "width": parsed_args["width"],
            #From pash log
            "exec_time": parsed_log["Execution time"],
            "backend_time": parsed_log["Backend time"],
            "compilation_time": parsed_log["Compilation time"],
            "preprocess_time": parsed_log["Preprocessing time"],
            "eager_nodes": parsed_log["Eager nodes"],
            "compiler_exit" : parsed_log["Compiler exited with code"],
            #From time
            "gnu_real": parsed_time["gnu_real"], 
            "gnu_usr": parsed_time["user"],
            "gnu_sys": parsed_time["sys"],
            "cpu%": parsed_time["cpu%"],
            "max_res_size": parsed_time["max_resident"],
            "avg_res_size": parsed_time["average_resident"],
            "major_pagefaults": parsed_time["major_pagefaults"],
            "minor_pagefaults": parsed_time["minor_pagefaults"],
        }

        #update local and global df
        df = df.append(data, ignore_index=True)
        self.df = self.df.append(data, ignore_index=True)

        return df
    def parse_file(self, file_path: str)->pd.DataFrame:
        """
        Parses a pa.sh log with path file_path
        Return:
            A single entry pandas dataframe
        """

    def parse_folder(self, path: str)->pd.DataFrame:
        """
        Parses all valid files ending with .log in the path directory.
        Params:
            path: should be a path to a directory containing .log files
        Return:
            pandas dataframe with all parsed logs
        """

    def get_df(self):
        return self.df

    def __parse_args__(self, args: str, args_of_interest) :
        lines = args.split("\n")
        args_dict = {i:False for i in args_of_interest}
        for line in lines:
            try:
                arg, val = line.split(" ")
                if arg in args_of_interest:
                    args_dict[arg] = val
            except:
                continue
                
        return args_dict
    
    def __parse_pash_log__(self, args: str, tags_of_interest) :
        lines = args.split("\n")
        log_dict = {i:False for i in tags_of_interest}
        
        for line in lines:
            try:
                tag, val = line.split(": ")
                if tag in tags_of_interest:
                    log_dict[tag] = val
            except:
                continue
        return log_dict

    def __parse_time_log__(self, timelog: str):
        data_start = timelog.find("Command being timed: ")
        time_data = timelog[data_start: ]

        lines = [line.split(": ")[1] for line in time_data.replace("\t", "").split("\n")[:-1]]
        if len(lines) < 23:
            lines = [False]*23
        data = {
            "command" : lines[0],
            "user" : lines[1],
            "sys" : lines[2],
            "cpu%" : lines[3],
            "gnu_real" : lines[4],
            "max_resident" : lines[9],
            "average_resident": lines[10],
            "major_pagefaults" : lines[11],
            "minor_pagefaults" : lines[12],
            "exit_status" : lines[22]
        }
        return data

#can be used in case we only can parse the time (default commands)
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