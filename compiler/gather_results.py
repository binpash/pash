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
UNIX50_RESULTS = "../evaluation/results/unix50/"

experiments = ["minimal_grep",
               "minimal_sort",
               "topn",
               "wf",
               "grep",
               "spell",
               "shortest_scripts",
               "diff"]

pretty_names = {"minimal_grep" : "grep",
                "minimal_sort" : "sort",
                "wf" : "wf",
                "topn" : "top-n",
                "grep" : "grep-old",
                "bigram" : "bi-gram",
                "spell" : "spell",
                "shortest_scripts" : "shortest-scripts",
                "diff" : "diff"}

def get_experiment_files(experiment, results_dir):
    files = [f for f in os.listdir(results_dir) if f.startswith(experiment)]
    return [int(f.split(experiment + "_")[1].split("_")[0]) for f in files]

def read_total_time(filename):
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

def read_distr_execution_time(filename):
    try:
        f = open(filename)
        times = []
        for line in f:
            if(line.startswith("Execution time")):
                milliseconds = line.split(": ")[1].split(" ")[0]
                times.append(float(milliseconds))
        f.close()
        return sum(times)
    except:
        print("!! WARNING: Filename:", filename, "not found!!!")
        return 0

def read_distr_total_compilation_time(filename):
    try:
        f = open(filename)
        times = []
        for line in f:
            if(line.startswith("Compilation time")
               or line.startswith("Optimization time")
               or line.startswith("Backend time")):
                milliseconds = line.split(": ")[1].split(" ")[0]
                times.append(float(milliseconds))
        f.close()
        # print(times)
        return sum(times)
    except:
        print("!! WARNING: Filename:", filename, "not found!!!")
        return 0

def collect_experiment_scaleup_times(prefix, scaleup_numbers):
    seq_numbers = [read_total_time('{}{}_seq.time'.format(prefix, n))
                   for n in scaleup_numbers]
    distr_numbers = [read_distr_execution_time('{}{}_distr.time'.format(prefix, n))
                     for n in scaleup_numbers]
    compile_numbers = [read_distr_total_compilation_time('{}{}_distr.time'.format(prefix, n))
                       for n in scaleup_numbers]
    # distr_numbers = [read_time('{}{}_distr.time'.format(prefix, n)) for n in all_scaleup_numbers]
    # cat_numbers = [read_time('{}{}_cat_distr.time'.format(prefix, n)) for n in all_scaleup_numbers]
    # compile_numbers = [read_time('{}{}_compile_distr.time'.format(prefix, n)) for n in all_scaleup_numbers]
    return (seq_numbers, distr_numbers, compile_numbers)

def collect_experiment_speedups(prefix, scaleup_numbers):
    seq_numbers, distr_numbers, compile_numbers = collect_experiment_scaleup_times(prefix, scaleup_numbers)
    # print(scaleup_numbers)
    # print(seq_numbers)
    # print(distr_numbers)
    # print(compile_numbers)
    distr_speedup = [seq_numbers[i] / t for i, t in enumerate(distr_numbers)]
    compile_distr_speedup = [seq_numbers[i] / (t + compile_numbers[i]) for i, t in enumerate(distr_numbers)]
    return (distr_speedup, compile_distr_speedup)

def collect_scaleup_times(experiment, results_dir):
    print(experiment)

    all_scaleup_numbers = list(set(get_experiment_files(experiment, results_dir)))
    all_scaleup_numbers.sort()
    all_scaleup_numbers = [i for i in all_scaleup_numbers if i > 1]
    prefix = '{}/{}_'.format(results_dir, experiment)
    distr_speedup, compile_distr_speedup = collect_experiment_speedups(prefix, all_scaleup_numbers)

    fig, ax = plt.subplots()

    ## Plot speedup
    ax.set_ylabel('Speedup')
    ax.set_xlabel('Level of Parallelism')
    # total_distr_speedup = [seq_numbers[i] / (t + compile_numbers[i] + cat_numbers[i]) for i, t in enumerate(distr_numbers)]
    # ax.plot(all_scaleup_numbers, seq_numbers, '-o', linewidth=0.5, label='Sequential')
    # ax.plot(all_scaleup_numbers, distr_numbers, '-o', linewidth=0.5, label='Distributed')
    ax.plot(all_scaleup_numbers, distr_speedup, '-o', linewidth=0.5, label='Parallel')
    ax.plot(all_scaleup_numbers, compile_distr_speedup, '-*', linewidth=0.5, label='+ Compile')
    # ax.plot(all_scaleup_numbers, total_distr_speedup, '-^', linewidth=0.5, label='+ Merge')
    # ax.plot(all_scaleup_numbers, all_scaleup_numbers, '-', color='tab:gray', linewidth=0.5, label='Ideal')


    # plt.yscale("log")
    plt.xticks(all_scaleup_numbers[1:])
    plt.legend(loc='lower right')
    plt.title(pretty_names[experiment])


    plt.tight_layout()
    plt.savefig(os.path.join('../evaluation/plots', "{}_throughput_scaleup.pdf".format(experiment)))

def collect_unix50_pipeline_scaleup_times(pipeline_number, unix50_results_dir, scaleup_numbers):
    prefix = '{}/unix50_pipeline_{}_'.format(unix50_results_dir, pipeline_number)
    return collect_experiment_speedups(prefix, scaleup_numbers)

def aggregate_unix50_results(all_results, scaleup_numbers):
    avg_distr_results = [[] for _ in scaleup_numbers]
    for pipeline in all_results:
        pipeline_distr_results = pipeline[0]
        # print(pipeline_distr_results)
        for i in range(len(scaleup_numbers)):
            avg_distr_results[i].append(pipeline_distr_results[i])

    for i in range(len(pipeline_distr_results)):
        avg_distr_results[i] = sum(avg_distr_results[i]) / len(avg_distr_results[i])

    return avg_distr_results

def collect_unix50_scaleup_times(unix50_results_dir):
    files = [f for f in os.listdir(unix50_results_dir)]
    # print(files)
    pipeline_numbers = sorted(list(set([f.split('_')[2] for f in files])))
    # print(pipeline_numbers)

    scaleup_numbers = [2, 4, 10, 20, 50]

    all_results = [collect_unix50_pipeline_scaleup_times(pipeline_number,
                                                         unix50_results_dir,
                                                         scaleup_numbers)
                   for pipeline_number in pipeline_numbers]
    # print(all_results)

    ## Plot individual speedups
    parallelism = 20
    individual_results = [distr_exec_speedup[scaleup_numbers.index(parallelism)]
                          for distr_exec_speedup, _ in all_results]
    print("Unix50 individual speedups for {} parallelism:".format(parallelism), individual_results)

    fig, ax = plt.subplots()
    ## Plot speedup
    ax.set_ylabel('Speedup')
    ax.set_xlabel('Pipeline')
    plt.bar(range(len(individual_results)), individual_results, align='center')
    plt.hlines([1], -1, len(individual_results) + 1, linewidth=0.8)
    plt.xlim(-1, len(individual_results) + 1)
    # plt.yscale("log")
    plt.yticks(range(1, 18, 2))
    # plt.ylim((0.1, 20))
    # plt.legend(loc='lower right')
    plt.title("Unix50 Individual Speedups")
    plt.tight_layout()
    plt.savefig(os.path.join('../evaluation/plots', "unix50_individual_speedups_{}.pdf".format(parallelism)))

    avg_results = aggregate_unix50_results(all_results, scaleup_numbers)
    print("Unix50 average speedup:", avg_results)

    ## Plot average speedup
    fig, ax = plt.subplots()

    ## Plot speedup
    ax.set_ylabel('Speedup')
    ax.set_xlabel('Level of Parallelism')
    ax.plot(scaleup_numbers, avg_results, '-o', linewidth=0.5, label='Parallel')
    # plt.yscale("log")
    plt.xticks(scaleup_numbers)
    plt.legend(loc='lower right')
    plt.title("Unix50 Throughput")
    plt.tight_layout()
    plt.savefig(os.path.join('../evaluation/plots', "unix50_throughput_scaleup.pdf"))


for experiment in experiments:
    collect_scaleup_times(experiment, RESULTS)

collect_unix50_scaleup_times(UNIX50_RESULTS)
