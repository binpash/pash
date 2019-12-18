import sys
import os
import numpy as np
import matplotlib.pyplot as plt

SMALL_SIZE = 16
MEDIUM_SIZE = 18
BIGGER_SIZE = 20

plt.rc('font', size=MEDIUM_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=BIGGER_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'


RESULTS = "../evaluation/results/"

experiments = ["minimal_grep",
               "minimal_sort",
               "topn",
               "wf",
               "grep",
               # "spell",
               "shortest_scripts"]

pretty_names = {"minimal_grep" : "grep",
                "minimal_sort" : "sort",
                "wf" : "wf",
                "topn" : "top-n",
                "grep" : "grep-old",
                "bigram" : "bi-gram",
                "spell" : "spell",
                "shortest_scripts" : "shortest-scripts"}

def get_experiment_files(experiment, results_dir):
    files = [f for f in os.listdir(results_dir) if f.startswith(experiment)]
    return [int(f.split(experiment + "_")[1].split("_")[0]) for f in files]

def read_time(filename):
    try:
        time = 0
        f = open(filename)
        for line in f:
            if(line.startswith("real")):
                minutes, seconds_milli = line.split("\t")[1].split("s")[0].split("m")
                seconds, milliseconds = seconds_milli.split(".")
                time = (int(minutes) * 60 + int(seconds)) * 1000 + int(milliseconds)
        f.close()
        return time
    except:
        print("!! WARNING: Filename:", filename, "not found!!!")
        return 0

def collect_scaleup_times(experiment, results_dir):
    print(experiment)

    all_scaleup_numbers = list(set(get_experiment_files(experiment, results_dir)))
    all_scaleup_numbers.sort()
    all_scaleup_numbers = [i for i in all_scaleup_numbers if i > 1]
    prefix = '{}/{}_'.format(results_dir, experiment)
    seq_numbers = [read_time('{}{}_seq.time'.format(prefix, n)) for n in all_scaleup_numbers]
    distr_numbers = [read_time('{}{}_distr.time'.format(prefix, n)) for n in all_scaleup_numbers]
    cat_numbers = [read_time('{}{}_cat_distr.time'.format(prefix, n)) for n in all_scaleup_numbers]
    compile_numbers = [read_time('{}{}_compile_distr.time'.format(prefix, n)) for n in all_scaleup_numbers]
    print(all_scaleup_numbers)
    print(seq_numbers)
    print(distr_numbers)
    fig, ax = plt.subplots()

    ## Plot speedup
    ax.set_ylabel('Speedup')
    ax.set_xlabel('Level of Parallelism')
    distr_speedup = [seq_numbers[i] / t for i, t in enumerate(distr_numbers)]
    compile_distr_speedup = [seq_numbers[i] / (t + compile_numbers[i]) for i, t in enumerate(distr_numbers)]
    total_distr_speedup = [seq_numbers[i] / (t + compile_numbers[i] + cat_numbers[i]) for i, t in enumerate(distr_numbers)]
    # ax.plot(all_scaleup_numbers, seq_numbers, '-o', linewidth=0.5, label='Sequential')
    # ax.plot(all_scaleup_numbers, distr_numbers, '-o', linewidth=0.5, label='Distributed')
    ax.plot(all_scaleup_numbers, distr_speedup, '-o', linewidth=0.5, label='Distributed')
    ax.plot(all_scaleup_numbers, compile_distr_speedup, '-*', linewidth=0.5, label='+ Compile')
    ax.plot(all_scaleup_numbers, total_distr_speedup, '-^', linewidth=0.5, label='+ Merge')
    # ax.plot(all_scaleup_numbers, all_scaleup_numbers, '-', color='tab:gray', linewidth=0.5, label='Ideal')


    # plt.yscale("log")
    plt.xticks(all_scaleup_numbers[1:])
    plt.legend(loc='lower right')
    plt.title(pretty_names[experiment])


    plt.tight_layout()
    plt.savefig(os.path.join('../evaluation/plots', "{}_throughput_scaleup.pdf".format(experiment)))

for experiment in experiments:
    collect_scaleup_times(experiment, RESULTS)
