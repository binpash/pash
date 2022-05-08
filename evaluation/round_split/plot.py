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

#This function only compares between folders
#dfs is tuple: (df_name, pd.df)
def plot_barchart_against(dfs: list, y_axis="gnu_real", test_names=None, save_file="output.png"):
    print(test_names)
    df1 = dfs[0][1]
    if test_names is not None:
        df1 = df1[df1["test_name"].isin(test_names)]
        labels = sorted(test_names)
    else:
        labels = sorted(df1.test_name.unique())

    x = np.arange(len(labels))  # the label locations
    width = 1/(len(dfs))*4/5  # the width of the bars

    rects = []
    for df_name, df in dfs:
        # Set filters here
        # df = df[(df.no_eager == False) & (df.dgsh_tee == True)]
        rects.append(df)

    bar_height = [[] for i in range(len(rects))]
    error_height = [[] for i in range(len(rects))]

    for label in labels:
        for i, rect in enumerate(rects):
            try:
                bar_height[i].append(rect.loc[rect["test_name"] == label][y_axis].mean())
                error_height[i].append(rect.loc[rect["test_name"] == label][y_axis].std())
            except:
                print("plot_shart_against: failed: ", label)
        

    fig, ax = plt.subplots(figsize=(20,10))
    for i in range(len(dfs)):
        label = ""
        df_name = dfs[i][0]
        offset = (i -  int(len(dfs)/2))
        if len(dfs)%2 == 0:
            offset += 1/2
        print(offset)
        ax.bar(x + offset*width, bar_height[i], width, yerr=error_height[i], label=label+df_name)

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
    log_dfs = []
    if len(argv) > 1:
        for folder in argv[1:]:
            log_parser = LogParser()
            df = log_parser.parse_folder(folder)
            log_dfs.append((folder, df))
    else:
        print("Usage: python3 plot.py [LOG FOLDER1] [LOG FOLDER2]...")
        os.exit(1)
    plot_barchart_against(log_dfs, test_names=['spell', 'set-diff', 'top-n', 'sort', 'bi-grams', 'sort-sort', 'diff', 'wf'], save_file="compare" + ".png")
