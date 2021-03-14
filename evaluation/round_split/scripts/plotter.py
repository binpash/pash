from log_parser import LogParser, DEFAULT_LOG_FOLDER
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def plot_barchart(df, column="exec_time"):
    # ax = df.plot.bar(x="test_name", y='exec_time', rot=0)
    # ax.savefig("barchart.png")
    # fig, ax = plt.subplots(figsize=(10,4))
    # for key, grp in df.groupby(['test_name']):
    #     ax.plot.bar(grp[column], label=key)
    # # fig = ax.get_figure()
    # fig.savefig('output.png')
    # ax.legend()
    # fi\
    # df = df[["test_name", "exec_time", "no_eager", "split_type"]]
    labels = df.test_name.unique()
    x = np.arange(len(labels))  # the label locations
    width = 0.2  # the width of the bars
    # rects2 = ax.bar(x + width/2, women_means, width, label='Women'
    
    r_split = df[(df.split_type == "r-split") & (df.no_eager == False)].set_index(labels)
    r_split_no_eager = df[(df.split_type == "r-split") & (df.no_eager == True)].set_index(labels)
    auto_split = df[(df.split_type == "auto-split") & (df.no_eager == False)].set_index(labels)
    auto_split_no_eager = df[(df.split_type == "auto-split") & (df.no_eager == True)].set_index(labels)
    
    r = []
    rn = []
    a = []
    an = []

    for label in labels:
        r.append(r_split.loc[r_split["test_name"] == label]["exec_time"][0])
        rn.append(r_split_no_eager.loc[r_split_no_eager["test_name"] == label]["exec_time"][0])
        a.append(auto_split.loc[auto_split["test_name"] == label]["exec_time"][0])
        an.append(auto_split_no_eager.loc[auto_split_no_eager["test_name"] == label]["exec_time"][0])
        
    # print(r_split)
    # print(r_split_no_eager)
    # print(auto_split)
    # print(auto_split_no_eager)
    fig, ax = plt.subplots(figsize=(20,10))
    rects1 = ax.bar(x - 3/2*width, r, width, label='r-split')
    rects2 = ax.bar(x - 1/2*width, rn, width, label='r-split-no-eager')
    rects3 = ax.bar(x + 1/2*width, a, width, label='auto-split')
    rects4 = ax.bar(x + 3/2*width, an, width, label='auto-split-no-eager')

    # rects2 = ax.bar(x + width/2, women_means, width, label='Women')

    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('ms')
    # ax.set_title('')
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    plt.xticks(rotation=90)
    ax.legend()
    # # ax = df.plot.bar()
    # ax = df.plot(kind='bar', figsize=(10,6), color="indigo", fontsize=13, x="test_name")
    # # ax.set_alpha(0.8)
    # fig = ax.get_figure()
    # # ax.get_figure().savefig("output.png")
    # # fig = ax.get_figure()
    fig.savefig('1G_width_4.png')

if __name__ == '__main__':
    #sample execution
    log_parser = LogParser()
    df = log_parser.parse_folder("tmp_1G_width4")
    plot_barchart(log_parser.get_df())