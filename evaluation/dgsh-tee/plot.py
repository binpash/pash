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
    print(test_names)
    if test_names is not None:
        df = df[df["test_name"].isin(test_names)]
        labels = sorted(test_names)
    else:
        labels = sorted(df.test_name.unique())

    x = np.arange(len(labels))  # the label locations
    width = 0.2  # the width of the bars
    
    r_split = df[(df.split_type == "r-split") & (df.no_eager == False) & (df.dgsh_tee == False)]

    r_split_dgsh_tee = df[(df.split_type == "r-split") & (df.no_eager == False) & (df.dgsh_tee == True)]

    auto_split = df[(df.split_type == "auto-split") & (df.no_eager == False) & (df.dgsh_tee == False)]

    auto_split_dgsh_tee = df[(df.split_type == "auto-split") & (df.no_eager == False) & (df.dgsh_tee == True)]
    
    rects = [r_split, r_split_dgsh_tee, auto_split, auto_split_dgsh_tee]
    bar_height = [[] for i in range(len(rects))]
    error_height = [[] for i in range(len(rects))]

    for label in labels:
        for i, rect in enumerate(rects):
            try:
                bar_height[i].append(rect.loc[rect["test_name"] == label][y_axis].mean())
                error_height[i].append(rect.loc[rect["test_name"] == label][y_axis].std())
            except:
                print(label)
        

    fig, ax = plt.subplots(figsize=(20,10))
    rects1 = ax.bar(x - 3/2*width, bar_height[0], width, yerr=error_height[0], label='r-split')
    rects2 = ax.bar(x - 1/2*width, bar_height[1], width, yerr=error_height[1], label='r-split-dgsh-tee')
    rects3 = ax.bar(x + 1/2*width, bar_height[2], width, yerr=error_height[2], label='auto-split')
    rects4 = ax.bar(x + 3/2*width, bar_height[3], width, yerr=error_height[3], label='auto-split-dgsh-tee')

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
    dirs_to_plot = argv[1:] if len(argv) > 1 else DEFAULT_LOG_FOLDER
    for d in dirs_to_plot:
        log_parser.parse_folder(d)
    # parse more directories if needed
    # df = log_parser.parse_folder("tmp_width4_10M")
    labels = ['minimal_sort', 'topn', 'wf', 'diff', 'set-diff', 'double_sort', 'sort']
    plot_barchart(log_parser.get_df(), test_names=labels, save_file=("_".join([os.path.dirname(d) for d in dirs_to_plot])) + ".png")
