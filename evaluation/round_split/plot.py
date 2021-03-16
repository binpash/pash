import sys
from test import PASH_TOP
sys.path.append(PASH_TOP)
from scripts.test_eval.logparser import LogParser, DEFAULT_LOG_FOLDER
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sys import argv
import os

def plot_barchart(df, y_axis="exec_time", test_names=None, save_file="output.png"):
    if test_names is not None:
        df = df[df["test_name"].isin(test_names)]
        labels = sorted(test_names)
    else:
        labels = sorted(df.test_name.unique())

    x = np.arange(len(labels))  # the label locations
    width = 0.2  # the width of the bars
    
    r_split = df[(df.split_type == "r-split") & (df.no_eager == False)]
    r_split_no_eager = df[(df.split_type == "r-split") & (df.no_eager == True)]
    auto_split = df[(df.split_type == "auto-split") & (df.no_eager == False)]
    auto_split_no_eager = df[(df.split_type == "auto-split") & (df.no_eager == True)]
    
    r = []
    rn = []
    a = []
    an = []
    error_r = []
    error_rn = []
    error_a = []
    error_an = []

    for label in labels:
        r.append(r_split.loc[r_split["test_name"] == label][y_axis].min())
        rn.append(r_split_no_eager.loc[r_split_no_eager["test_name"] == label][y_axis].min())
        a.append(auto_split.loc[auto_split["test_name"] == label][y_axis].min())
        an.append(auto_split_no_eager.loc[auto_split_no_eager["test_name"] == label][y_axis].min())

        error_r.append(r_split.loc[r_split["test_name"] == label][y_axis].max() - r[-1])
        error_rn.append(r_split_no_eager.loc[r_split_no_eager["test_name"] == label][y_axis].max() - rn[-1])
        error_a.append(auto_split.loc[auto_split["test_name"] == label][y_axis].max() - a[-1])
        error_an.append(auto_split_no_eager.loc[auto_split_no_eager["test_name"] == label][y_axis].max() - an[-1])
        

    fig, ax = plt.subplots(figsize=(20,10))
    rects1 = ax.bar(x - 3/2*width, r, width, yerr=error_r, label='r-split')
    rects2 = ax.bar(x - 1/2*width, rn, width, yerr=error_rn, label='r-split-no-eager')
    rects3 = ax.bar(x + 1/2*width, a, width, yerr=error_a, label='auto-split')
    rects4 = ax.bar(x + 3/2*width, an, width, yerr=error_an, label='auto-split-no-eager')

    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('ms')
    # ax.set_title('')
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    plt.xticks(rotation=90)
    ax.legend()
    
    fig.savefig(save_file)

if __name__ == '__main__':
    #sample execution
    log_parser = LogParser()
    dir_to_plot = argv[1] if len(argv) > 1 else DEFAULT_LOG_FOLDER
    df = log_parser.parse_folder(dir_to_plot)
    # parse more directories if needed
    # df = log_parser.parse_folder("tmp2_1G_width4bbackup")

    plot_barchart(log_parser.get_df(), save_file=os.path.dirname(dir_to_plot) + ".png")