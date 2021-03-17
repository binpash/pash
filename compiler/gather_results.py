import copy
import sys
import math
import argparse
import os
import re
import statistics
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as pltlines
import matplotlib.ticker as plticker

parser = argparse.ArgumentParser(description='Produce plots from various experiments with PaSh.')
parser.add_argument('--eurosys2021',
                    action='store_true',
                    help='generates the plots for all the experiments in the EuroSys2021 paper')
parser.add_argument('--all',
                    action='store_true',
                    help='generates all plots')

args = parser.parse_args()

if args.all is True:
    args.eurosys2021 = True

if not args.all and not args.eurosys2021:
    print("You have to specify some plot to generate!")
    print("See command usage with --help.")
    exit(0)

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


MICROBENCHMARKS = "../evaluation/microbenchmarks/"
RESULTS = "../evaluation/results/"
SMALL_RESULTS = "../evaluation/results/eurosys_small/"
UNIX50_RESULTS = "../evaluation/results/unix50/"
SMALL_UNIX50_RESULTS = "../evaluation/results/unix50_4_1073741824/"
BIG_UNIX50_RESULTS = "../evaluation/results/unix50_16_10737418240/"
COARSE_UNIX50_RESULTS = "../evaluation/results/unix50-naive/"
OASA_RESULTS = "../evaluation/buses/results/"


all_experiments = ["minimal_grep",
                   "minimal_sort",
                   "topn",
                   "wf",
                   "grep",
                   "spell",
                   "shortest_scripts",
                   "diff",
                   "bigrams",
                   "alt_bigrams",
                   "set-diff",
                   "double_sort"]

pretty_names = {"minimal_grep" : "nfa-regex",
                "minimal_sort" : "sort",
                "wf" : "wf",
                "topn" : "top-n",
                "grep" : "filter",
                "bigrams" : "bi-grams",
                "alt_bigrams" : "bi-grams-opt",
                "spell" : "spell",
                "shortest_scripts" : "shortest-scripts",
                "diff" : "difference",
                "set-diff" : "set-difference",
                "double_sort" : "sort-sort"}

structures = {"minimal_grep" : "$3\\times\\tsta$",
              "minimal_sort" : "$\\tsta, \\tpur$",
              "wf" : "$3\\times\\tsta, 3\\times\\tpur$",
              "topn" : "$2\\times\\tsta, 4\\times\\tpur$",
              "grep" : "$3\\times\\tsta$",
              "bigrams" : "$3\\times\\tsta, 3\\times\\tpur$",
              "alt_bigrams" : "$3\\times\\tsta, \\tpur$",
              "spell" : "$4\\times\\tsta, 3\\times\\tpur$",
              "shortest_scripts" : "$5\\times\\tsta, 2\\times\\tpur$",
              "diff" : "$2\\times\\tsta, 3\\times\\tpur$",
              "set-diff" : "$5\\times\\tsta, 2\\times\\tpur$",
              "double_sort" : "$\\tsta, 2\\times\\tpur$"}

highlights = {"minimal_grep" : "complex NFA regex",
              "minimal_sort" : "\\tti{sort}ing",
              "wf" : "double \\tti{sort}, \\tti{uniq} reduction",
              "topn" : "double \\tti{sort}, \\tti{uniq} reduction",
              "grep" : "IO-intensive, computation-light",
              "bigrams" : "stream shifting and merging",
              "alt_bigrams" : "optimized version of bigrams",
              "spell" : "comparisons (\\tti{comm})",
              "shortest_scripts" : "long \\tsta pipeline ending with \\tpur",
              "diff" : "non-parallelizable \\tti{diff}ing",
              "set-diff" : "two pipelines merging to a \\tti{comm}",
              "double_sort" : "parallelizable \\tpur after \\tpur"}

input_filename_sizes = {"1G": "1~GB",
                        "3G": "3~GB",
                        "10G": "10~GB",
                        "100G": "100~GB",
                        "1M": "1~MB",
                        "10M": "10~MB",
                        "100M": "100~MB",
                        "all_cmds_x1000": "85~MB"}

suffix_to_runtime_config = {"distr": "eager",
                            "distr_auto_split": "split",
                            "distr_no_task_par_eager": "blocking-eager",
                            "distr_no_eager": "no-eager",
                            "distr_auto_split_fan_in_fan_out": "no-aux-cat-split"}

all_line_plots = ["split",
                  "mini-split",
                  "eager",
                  "blocking-eager",
                  "no-eager",
                  "no-aux-cat-split"]

file_suffixes = {"split": "distr_auto_split.time",
                 "mini-split": "distr_auto_split.time",
                 "eager": "distr.time",
                 "blocking-eager": "distr_no_task_par_eager.time",
                 "no-eager": "distr_no_eager.time",
                 "no-aux-cat-split": "distr_auto_split_fan_in_fan_out.time"}


class LinePlotConfig:
    ## TODO: Have a default label
    def __init__(self, linestyle, color, label, linewidth=0.5):
        self.linestyle = linestyle
        self.color = color
        self.label = label
        self.linewidth = linewidth
    
    def plot(self, xvalues, yvalues, ax):
        return ax.plot(xvalues, yvalues, self.linestyle, linewidth=self.linewidth, color=self.color, label=self.label)

    ## The following two methods are needed since the linestyle argument in plot
    ## does not correspond to the one in lines2D
    def get_ls(self):
        linestyle = self.linestyle
        if(len(linestyle) == 2 and linestyle[0] == '-' and not linestyle[1] == '-'):
            return '-'
        return linestyle

    def get_marker(self):
        linestyle = self.linestyle
        if(len(linestyle) == 2 and linestyle[0] == '-' and not linestyle[1] == '-'):
            return linestyle[1]
        return None

    def get_artist(self):
        marker = self.get_marker()
        linestyle = self.get_ls()
        return pltlines.Line2D([], [], marker=marker, linestyle=linestyle, linewidth=self.linewidth, color=self.color, label=self.label)

default_line_plot_configs = {'eager': LinePlotConfig('-D', 'tab:red', 'Parallel w/o split'),
                             'split': LinePlotConfig('-o', 'tab:blue', 'Parallel'),
                             'mini-split': LinePlotConfig('-o', 'tab:blue', 'Parallel'),
                             'blocking-eager': LinePlotConfig('-p', 'orange', 'Blocking Eager'),
                             'no-eager': LinePlotConfig('-^', 'green', 'No Eager'),
                             'no-aux-cat-split': LinePlotConfig('-v', 'brown', 'No Aux Cat-Split')}

class Config:
    def __init__(self, pash, width=None, runtime=None):
        if(pash):
            assert(Config.is_runtime_valid(runtime))
            self.pash = True
            self.width = width
            self.runtime = runtime
        else:
            self.pash = False
    
    def __repr__(self):
        if(self.pash):
            return("PaSh Config(width={}, runtime={})".format(self.width, self.runtime))
        else:
            return("No PaSh")

    def is_runtime_valid(runtime):
        possible_runtimes = ["no-eager",
                             "blocking-eager",
                             "eager",
                             "mini-split",
                             "split",
                             "no-aux-cat-split"]
        return runtime in possible_runtimes


class Result:
    def __init__(self, script, config, value, description=""):
        self.script = script
        self.config = config
        self.value = value
        self.description = description
    
    def __repr__(self):
        return("Result(description={}, script={}, value={}, config={})".format(self.description,
                                                                               self.script,
                                                                               self.value,
                                                                               self.config))

    def __truediv__(self, other):
        if not isinstance(other, (int, float, Result)):
            return NotImplemented

        ## TODO: Change that to Result too
        return self.value / other

    def __rtruediv__(self, other):
        if not isinstance(other, (int, float, Result)):
            return NotImplemented

        if(self.value == 0):
            print("Division by zero")
            return 0

        ## TODO: Change that to Result too
        return other / self.value

    def __add__(self, other):
        if not isinstance(other, (int, float, Result)):
            return NotImplemented

        if not self.description == "execution time":
            return NotImplemented

        ## TODO: Change that to Result too
        return other + self.value

    def __radd__(self, other):
        return self + other

class ResultVector:
    def __init__(self, script, config, results, description, xvalues, xaxis):
        self.script = script
        self.config = config
        self.results = results 
        self.description = description
        self.xvalues = xvalues
        self.xaxis = xaxis

    def __repr__(self):
        return("Results(description={}, script={}, config={})".format(self.description,
                                                                      self.script,
                                                                      self.value,
                                                                      self.config))

    def plot(self, ax, linestyle='--', linewidth=0.5, color=None, label=None):
        if(label is None):
            ## TODO: If a label is not passed, we can generate one through config
            return NotImplemented
        
        if(color is None):
            return ax.plot(self.xvalues, self.results, linestyle, linewidth=linewidth, label=label)
        else:
            ## This might not be a necessary case-split
            return ax.plot(self.xvalues, self.results, linestyle, linewidth=linewidth, color=color, label=label)

    def __iter__(self):
        ''' Returns the iterator for the list of results. '''
        return self.results.__iter__()

def safe_zero_div(a, b):
    if(a is None or b is None):
        return None
    elif(b == 0):
        print("WARNING: Division by zero")
        return 0
    else:
        return a / b


##
## Read From Files
##

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
        if(len(times) == 0):
            return 0
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

def read_distr_command_number(filename):
    try:
        with open(filename) as f:
            for line in f:
                if(line.startswith("Total nodes after optimization:")):
                    number_of_commands = int(line.split(": ")[1])
                    return number_of_commands
            return "\\todo{X}"
    except:
        print("!! WARNING: Filename:", filename, "not found!!!")
        return "\\todo{X}"

def check_output_diff_correctness_for_experiment(filename):
    try:
        with open(filename) as f:
            for line in f:
                if(line.startswith("Files ")
                   and line.rstrip().endswith("are identical")):
                    return "Correct"
            return "Wrong"
    except:
        # print("!! WARNING: Filename:", filename, "not found!!!")
        return "Not Found"

##
## Result Collection
##

def script_name_from_prefix(prefix):
    ## TODO: Move this outside
    script_name = prefix.split("/")[-1].rstrip("_")
    return script_name

def runtime_config_from_suffix(suffix):
    runtime = suffix_to_runtime_config[suffix.split(".")[0]]
    return runtime

def sequential_experiment_exec_time(prefix, scaleup_number):
    config = Config(pash=False)
    value = read_total_time('{}{}_seq.time'.format(prefix, scaleup_number))
    description = "execution time"
    script_name = script_name_from_prefix(prefix)
    result = Result(script_name, config, value, description)
    return result

def distributed_experiment_exec_time(prefix, scaleup_number, suffix):
    pash_runtime = runtime_config_from_suffix(suffix)
    # print(pash_runtime)
    config = Config(pash=True, width=scaleup_number, runtime=pash_runtime)
    value = read_distr_execution_time('{}{}_{}'.format(prefix, scaleup_number, suffix))
    description = "execution time"
    script_name = script_name_from_prefix(prefix)
    result = Result(script_name, config, value, description)
    return result

def collect_distr_experiment_execution_times(prefix, suffix, scaleup_numbers):
    numbers = [distributed_experiment_exec_time(prefix, n, suffix)
               for n in scaleup_numbers]

    ## TODO: Turn to Result
    compile_numbers = [read_distr_total_compilation_time('{}{}_{}'.format(prefix, n, suffix))
                       for n in scaleup_numbers]
    return (numbers, compile_numbers)

def collect_experiment_scaleup_times(prefix, scaleup_numbers, suffix="distr.time"):
    ## Since we have the same input size in all cases, only use the
    ## one sequential execution for the sequential time
    seq_number = sequential_experiment_exec_time(prefix, scaleup_numbers[0])
    distr_numbers, compile_numbers = collect_distr_experiment_execution_times(prefix, suffix, scaleup_numbers)
    return (seq_number, distr_numbers, compile_numbers)

def collect_baseline_experiment_speedups(prefix, scaleup_numbers, base_seq):
    seq_numbers = [sequential_experiment_exec_time(prefix, n)
                   for n in scaleup_numbers]
    speedup = [safe_zero_div(base_seq, t) for t in seq_numbers]
    return speedup

def collect_distr_experiment_speedup_seq_time(prefix, scaleup_numbers, suffix="distr.time"):
    distr_speedup, _, seq_time = collect_distr_experiment_speedup_with_compilation(prefix, scaleup_numbers, suffix)
    return (distr_speedup, seq_time)

def collect_distr_experiment_speedup(prefix, scaleup_numbers, suffix="distr.time"):
    distr_speedup, _ = collect_distr_experiment_speedup_seq_time(prefix, scaleup_numbers, suffix=suffix)
    return distr_speedup

def collect_distr_experiment_speedup_with_compilation(prefix, scaleup_numbers, suffix="distr.time"):
    seq_number, distr_numbers, compile_numbers = collect_experiment_scaleup_times(prefix, scaleup_numbers, suffix=suffix)
    distr_speedups = [safe_zero_div(seq_number, t) for i, t in enumerate(distr_numbers)]
    compile_distr_speedups = [safe_zero_div(seq_number, t + compile_numbers[i]) for i, t in enumerate(distr_numbers)]

    ## TODO: Populate the Results using a method in Results   
    ## TODO: Pass these as arguments to the other methods 
    ##       (or if populating through a method in Results we can just fill there) 
    # script_name = script_name_from_prefix(prefix)
    # description = "speedups"
    # runtime = runtime_config_from_suffix(suffix)
    # config = Config(pash=True, width="Var", runtime=runtime)
    # result_vec = ResultVector(script_name, config, distr_speedups, 
    #                           description, scaleup_numbers, "width")
    # return (result_vec, compile_distr_speedups)
    return (distr_speedups, compile_distr_speedups, seq_number)


def collect_experiment_command_number(prefix, suffix, scaleup_numbers):
    command_numbers = [read_distr_command_number('{}{}_{}'.format(prefix, n, suffix))
                       for n in scaleup_numbers]
    return command_numbers

def check_output_diff_correctness(prefix, scaleup_numbers):
    global all_line_plots
    global file_suffixes
    
    result_correctness = {}
    for line_plot in all_line_plots:
        suffix = file_suffixes[line_plot]
        result_correctness[line_plot] = [(n, check_output_diff_correctness_for_experiment('{}{}_{}'.format(prefix, n, suffix)))
                                         for n in scaleup_numbers]
    return result_correctness

def collect_scaleup_line_speedups(experiment, all_scaleup_numbers, results_dir):
    global all_line_plots
    global file_suffixes

    print("Collecting results for:", experiment)

    prefix = '{}/{}_'.format(results_dir, experiment)

    ## Gather results
    all_speedup_results = {}
    sequential_time = 0
    for line_plot in all_line_plots:
        file_suffix = file_suffixes[line_plot]
        try:
            speedup_results, sequential_time = collect_distr_experiment_speedup_seq_time(prefix,
                                                                                         all_scaleup_numbers,
                                                                                         file_suffix)
            all_speedup_results[line_plot] = speedup_results
        except ValueError:
            ## TODO: Should we do anything here?
            pass

    ## Check if outputs are correct
    output_diff = check_output_diff_correctness(prefix, all_scaleup_numbers)
    
    return (all_speedup_results, output_diff, sequential_time)


def plot_scaleup_lines(experiment, all_scaleup_numbers, all_speedup_results, custom_scaleup_plots, 
                       ax, line_plot_configs=default_line_plot_configs):
    print("Plotting:", experiment)

    default_line_plots = ["split",
                          "eager",
                          "blocking-eager",
                          "no-eager",
                          "no-aux-cat-split"]

    ## Decide which lines to plot
    line_plots = default_line_plots
    if(experiment in custom_scaleup_plots):
        line_plots = custom_scaleup_plots[experiment]

    ## Compute the best speedups (for averages)
    best_result = []
    if("no-eager" in line_plots):
        best_result = all_speedup_results["no-eager"]
    if("blocking-eager" in line_plots):
        best_result = all_speedup_results["blocking-eager"]
    if("eager" in line_plots):
        best_result = all_speedup_results["eager"]
    if("split" in line_plots):
        split_res = all_speedup_results["split"]
        if(sum(split_res) > sum(best_result)):
            best_result = all_speedup_results["split"]
    elif("mini-split" in line_plots):
        split_res = all_speedup_results["mini-split"]
        if(sum(split_res) > sum(best_result)):
            best_result = all_speedup_results["mini-split"]

    ## We need to return the no-eager speedups as the baseline non runtime primitives.
    no_eager_distr_speedup = all_speedup_results["eager"]
    if("no-eager" in line_plots):
        no_eager_distr_speedup = all_speedup_results["no-eager"]


    ## Compute if a split line exists to change the color of the non-split top line
    split_exists = "split" in line_plots

    lines = []
    to_plot_lines = [line_plot for line_plot in line_plots
                     if line_plot in all_speedup_results]
    maximum_y = 0
    for line_plot in to_plot_lines:
        ## Get the result
        speedup_results = all_speedup_results[line_plot]

        ## Aggregate all the speedups to find the maximum one
        maximum_y = max(maximum_y, max(speedup_results))

        ## Get the config
        plot_config = copy.deepcopy(line_plot_configs[line_plot])

        ## If a split exists, the eager should be red
        ## TODO: Move this logic outside. It should not be part of the plotting function
        if(line_plot == "eager" and not split_exists):
            plot_config.color = 'tab:blue'
            plot_config.linestyle = '-o'
            plot_config.label = 'Parallel'

        line, = plot_config.plot(all_scaleup_numbers, speedup_results, ax)
        lines.append(line)

    ## Set the ylim
    ## TODO: Move this outside
    ax.set_ylim(top=maximum_y*1.15)

    return lines, best_result, no_eager_distr_speedup

def plot_sort_with_baseline(results_dir):

    all_scaleup_numbers = [2, 4, 8, 16, 32, 64]
    sort_prefix = '{}/sort_'.format(results_dir)
    baseline_sort_prefix = '{}/baseline_sort/baseline_sort_'.format(results_dir)
    baseline_sort_opt_prefix = '{}/baseline_sort/baseline_sort_opt_'.format(results_dir)

    ## Collect all sort numbers
    seq_number, distr_numbers, _ = collect_experiment_scaleup_times(sort_prefix, all_scaleup_numbers)
    sort_distr_speedup = [safe_zero_div(seq_number, t) for i, t in enumerate(distr_numbers)]
    # sort_distr_speedup = collect_distr_experiment_speedup(sort_prefix, all_scaleup_numbers)
    baseline_sort_distr_speedup = collect_baseline_experiment_speedups(baseline_sort_prefix,
                                                                       [1] + [num*2
                                                                              for num in all_scaleup_numbers],
                                                                       seq_number)
    # baseline_sort_opt_distr_speedup = collect_baseline_experiment_speedups(baseline_sort_opt_prefix,
    #                                                                        [1] + all_scaleup_numbers[1:],
    #                                                                        seq_numbers[0])

    # output_diff = check_output_diff_correctness(prefix, all_scaleup_numbers)

    fig, ax = plt.subplots()

    ## Plot speedup
    ax.set_ylabel('Speedup')
    ax.set_xlabel('--width')
    ax.plot(all_scaleup_numbers, sort_distr_speedup, '-o', linewidth=0.5, label='Pash')
    ## Add the no eager times if they exist
    # try:
    #     no_task_par_eager_distr_speedup = collect_distr_experiment_speedup(sort_prefix,
    #                                                                         all_scaleup_numbers,
    #                                                                         'distr_no_task_par_eager.time')
    #     ax.plot(all_scaleup_numbers, no_task_par_eager_distr_speedup, '-p', linewidth=0.5, label='Pash - Blocking Eager')
    # except ValueError:
    #     pass

    try:
        no_eager_distr_speedup = collect_distr_experiment_speedup(sort_prefix,
                                                                  all_scaleup_numbers,
                                                                  'distr_no_eager.time')
        ax.plot(all_scaleup_numbers, no_eager_distr_speedup, '-^', linewidth=0.5, label='Pash - No Eager')
    except ValueError:
        pass

    ax.plot(all_scaleup_numbers, baseline_sort_distr_speedup[1:], '-p', linewidth=0.5, label='sort --parallel')
    # ax.plot(all_scaleup_numbers, baseline_sort_opt_distr_speedup[1:], '-', linewidth=0.5, label='sort --parallel -S 30%')


    plt.xticks(all_scaleup_numbers[1:])
    plt.legend(loc='lower right')
    # plt.title("Comparison with sort --parallel")


    plt.tight_layout()
    plt.savefig(os.path.join('../evaluation/plots', "sort_baseline_comparison_scaleup.pdf"),bbox_inches='tight')


def collect_format_input_size(experiment):
    raw_size = collect_input_size(experiment)
    try:
        result = input_filename_sizes[raw_size]
    except:
        result = "\\todo{UNKNOWN}"
    return result

def collect_input_size(experiment):
    env_file = os.path.join(MICROBENCHMARKS, '{}_env.sh'.format(experiment))
    with open(env_file) as file:
        input_file_names = [line.rstrip().split("=")[1] for line in file.readlines() if line.startswith("IN")]
        assert(len(input_file_names) == 1)
        input_file_name = input_file_names[0]
    # print(input_file_name)
    # try:
    #     input_size = os.stat(input_file_name).st_size
    # except:
    #     input_size = 0
    clean_name = input_file_name.split('/')[-1].split('.')[0]
    return clean_name

def separate_into_lines(script):
    indent_len = 4
    lines = [script]
    while (len(lines[-1]) > 90):
        last_line = lines[-1]
        last_break_index = -1
        for i in range(90):
            if(last_line[i] in ["&", "|"]):
                last_break_index = i

        if(last_break_index == -1):
            break
        new_pre_last_line = last_line[:(last_break_index+1)]
        new_last_line = " " * indent_len + last_line[(last_break_index+1):].lstrip()
        lines[-1] = new_pre_last_line
        lines.append(new_last_line)
    return lines


def collect_script(experiment):
    # print(experiment)
    script_file = os.path.join(MICROBENCHMARKS, '{}.sh'.format(experiment))
    with open(script_file) as file:
        script_lines = file.readlines()
    clean_script_lines = [line.lstrip() for line in script_lines
                            if not line.lstrip().startswith('#')]
    clean_script_lines = [line.split('#')[0].rstrip(' ') for line in clean_script_lines]
    ## Remove mkfifo, cat, rm
    clean_script_lines = [line for line in clean_script_lines
                            if not line.startswith("mkfifo") and
                                not line.startswith("rm")]
    script = "".join(clean_script_lines).rstrip()
    script = script.replace("|\n", "| ")
    script = script.replace("&\n", "& ")
    script = script.replace("\n", "; ")
    print(len(script), script)
    lines = separate_into_lines(script)
    wrapped_lines = ['\lstinline[columns=fixed,basicstyle=\\footnotesize\\ttfamily]!{}!'.format(line)
                     for line in lines]
    wrapped_script = '\\\\ ~&~ '.join(wrapped_lines)
    return wrapped_script

def format_time_seconds(time_milliseconds):
    time_seconds = time_milliseconds / 1000
    time_minutes = int(time_seconds // 60)
    time_only_seconds = time_seconds % 60
    if(time_minutes > 0):
        formatted_time = '{}m{:.3f}s'.format(time_minutes, time_only_seconds)
    else:
        formatted_time = '{:.3f}s'.format(time_seconds)
    return formatted_time

def format_time_milliseconds(time_milliseconds):
    time_seconds = int(time_milliseconds // 1000)
    time_only_milliseconds = time_milliseconds % 1000
    if(time_seconds > 0):
        formatted_time = '{}s{:.3f}ms'.format(time_seconds, time_only_milliseconds)
    else:
        formatted_time = '{:.3f}ms'.format(time_milliseconds)
    return formatted_time


def generate_table_header(full=True):
    header = []
    if(not full):
        # header += ['\\begin{tabular*}{\\textwidth}{l @{\\extracolsep{\\fill}} llll}']
        header += ['\\begin{tabular*}{\\textwidth}{l @{\\extracolsep{\\fill}} l}']
        header += ['\\toprule']
        # header += ['Script ~&~ Structure & Input &'
        #            'Seq. Time & Highlights \\\\']
        header += ['Name ~&~ Script \\\\']
    else:
        header += ['\\begin{tabular*}{\\textwidth}{l @{\\extracolsep{\\fill}} llllllll}']
        header += ['\\toprule']
        header += ['Script ~&~ Structure & Input &'
                   'Seq. Time & \\multicolumn{2}{l}{\\#Nodes(16, 64)} &'
                   '\\multicolumn{2}{l}{Compile Time (16, 64)} & Highlights \\\\']
    header += ['\\midrule']
    return "\n".join(header)

def generate_table_footer(full=True):
    footer = []
    footer += ['\\bottomrule']
    footer += ['\\end{tabular*}']
    return "\n".join(footer)


def generate_experiment_line(experiment, results_dir=RESULTS, full=True, small=False):
    line = []
    line += [pretty_names[experiment], '~&~']

    if(not full):
        script = collect_script(experiment)
        line += [script, '\\\\']
        return " ".join(line)

    line += [structures[experiment], '&']
    ## Collect and output the input size
    input_size = collect_format_input_size(experiment)
    line += [input_size, '&']

    if(small is False):
        suffix='distr.time'
        if(experiment in ["spell", "bigrams", "double_sort"]):
            suffix='distr_auto_split.time'
    else:
        suffix='distr_auto_split.time'

    ## Collect and output the sequential time for the experiment
    scaleup_numbers = [2, 16, 64]
    experiment_results_prefix = '{}/{}_'.format(results_dir, experiment)
    seq_time, _, compile_times = collect_experiment_scaleup_times(experiment_results_prefix, scaleup_numbers, suffix=suffix)
    seq_time_seconds = format_time_seconds(seq_time)
    # seq_time_seconds = seq_times[0] / 1000
    line += [seq_time_seconds, '&']

    commands_16, commands_64 = collect_experiment_command_number(experiment_results_prefix,
                                                                suffix, [16, 64])
    if(full):
        line += ['{} & {}'.format(commands_16, commands_64), '&']

    ## Collect and output compile times
    compile_time_16_milliseconds = compile_times[1]
    compile_time_64_milliseconds = compile_times[2]
    if(full):
        line += ['{} & {}'.format(format_time_seconds(compile_time_16_milliseconds),
                                  format_time_seconds(compile_time_64_milliseconds)), '&']
    line += [highlights[experiment], '\\\\']
    return " ".join(line)

def generate_tables(experiments, results_dir=RESULTS, table_suffix="", small=False):
    generate_tex_table(experiments, results_dir=results_dir, table_suffix=table_suffix, small=small)

def generate_tex_table(experiments, results_dir=RESULTS, table_suffix="", small=False):
    header = generate_table_header()
    lines = []
    for experiment in experiments:
        line = generate_experiment_line(experiment, results_dir=results_dir, small=small)
        # print(line)
        lines.append(line)
    data = "\n".join(lines)
    footer = generate_table_footer()
    table_tex = "\n".join([header, data, footer])
    tex_filename = os.path.join('../evaluation/plots', 'microbenchmarks-table{}.tex'.format(table_suffix))
    with open(tex_filename, 'w') as file:
        file.write(table_tex)

def generate_tex_coarse_table(experiments):
    header = generate_table_header(full=False)
    lines = []
    for experiment in experiments:
        line = generate_experiment_line(experiment, full=False)
        # print(line)
        lines.append(line)
    data = "\n".join(lines)
    footer = generate_table_footer(full=False)
    table_tex = "\n".join([header, data, footer])
    tex_filename = os.path.join('../evaluation/plots', 'coarse-microbenchmarks-table.tex')
    with open(tex_filename, 'w') as file:
        file.write(table_tex)


def collect_unix50_pipeline_scaleup_times(pipeline_number, unix50_results_dir, scaleup_numbers, suffix='distr.time'):
    prefix = '{}/unix50_pipeline_{}_'.format(unix50_results_dir, pipeline_number)
    distr_speedups = collect_distr_experiment_speedup(prefix, scaleup_numbers, suffix)
    absolute_seq_time, _, _ = collect_experiment_scaleup_times(prefix, scaleup_numbers)
    return (distr_speedups, absolute_seq_time)

def collect_unix50_pipeline_coarse_scaleup_times(pipeline_number, unix50_results_dir, scaleup_numbers):
    prefix = '{}/unix50_pipeline_{}_'.format(unix50_results_dir, pipeline_number)
    # distr_speedups = collect_distr_experiment_speedup(prefix, scaleup_numbers)
    no_eager_distr_speedup = collect_distr_experiment_speedup(prefix,
                                                              scaleup_numbers,
                                                              'distr_no_eager.time')
    absolute_seq_times = [sequential_experiment_exec_time(prefix, scaleup_numbers[0])
                          for _ in scaleup_numbers]
    return (no_eager_distr_speedup, absolute_seq_times)

def collect_unix50_pipeline_fan_in_fan_out_scaleup_times(pipeline_number, unix50_results_dir, scaleup_numbers):
    prefix = '{}/unix50_pipeline_{}_'.format(unix50_results_dir, pipeline_number)
    fan_in_fan_out_distr_speedup = collect_distr_experiment_speedup(prefix,
                                                                    scaleup_numbers,
                                                                    'distr_auto_split_fan_in_fan_out.time')
    absolute_seq_times = [sequential_experiment_exec_time(prefix, scaleup_numbers[0])
                          for _ in scaleup_numbers]
    return (fan_in_fan_out_distr_speedup, absolute_seq_times)


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

def compute_and_print_aggrs(individual_results, absolute_seq_times_s):
    mean = sum(individual_results) / len(individual_results)
    median = statistics.median(individual_results)
    geo_mean = math.exp(np.log(individual_results).sum() / len(individual_results))
    weighted_res = [i*a for i, a in zip(individual_results, absolute_seq_times_s)]
    weighted_avg = sum(weighted_res) / sum(absolute_seq_times_s)
    print("  Mean:", mean)
    print("  Median:", median)
    print("  Geometric Mean:", geo_mean)
    print("  Weighted Average:", weighted_avg)
    return mean, median, geo_mean, weighted_res, weighted_avg

def make_unix50_bar_chart(all_results, scaleup_numbers, parallelism, small_prefix=""):

    ## Sort by speedup
    sorted_indices = sorted(range(len(all_results)), key=lambda x:all_results[x][0][scaleup_numbers.index(parallelism)], reverse=True)
    sorted_all_results = sorted(all_results, key=lambda x: x[0][scaleup_numbers.index(parallelism)], reverse=True)
    print("Unix50 Par: {} Sorted Indices:".format(parallelism))
    print("|------------------------") 
    indices_map = sorted([(v,i) for i, v in enumerate(sorted_indices)], key=lambda x: x[0])
    print("\n".join(["Old: {}, New: {}".format(v,i) for v, i in indices_map]))
    print("|------------------------") 
    ## Filter small exec times.
    sorted_all_results = [res for res in sorted_all_results
                          if (res[1] / 1000) > 0.1]

    ## Plot individual speedups
    individual_results = [distr_exec_speedup[scaleup_numbers.index(parallelism)]
                          for distr_exec_speedup, _ in sorted_all_results]
    absolute_seq_times_s = [absolute_seq / 1000
                            for _, absolute_seq in sorted_all_results]
    print("Unix50 individual speedups for {} parallelism:".format(parallelism), individual_results)

    aggrs = compute_and_print_aggrs(individual_results, absolute_seq_times_s)
    mean = aggrs[0]

    w = 0.2
    ind = np.arange(len(individual_results))
    speedup_color = 'tab:blue'
    seq_time_color = 'tab:red'



    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
    # print(fig.get_size_inches())
    fig.set_size_inches(14, 6)

    ## Plot line at 1
    plt.xlim(-1, len(individual_results))

    ## Plot speedup
    ax2.set_ylabel('Speedup', color=speedup_color)
    ax2.bar(ind-w, individual_results, width=2*w, align='center', color=speedup_color)
    ax2.tick_params(axis='y', labelcolor=speedup_color)

    ax1.set_ylabel('Sequential Time (s)', color=seq_time_color)
    ax1.set_yscale("log")
    ax1.set_ylim([5, 6000])
    ax1.set_xlabel('Script Index')
    ax1.bar(ind+w, absolute_seq_times_s, width=2*w, align='center', color=seq_time_color)
    ax1.tick_params(axis='y', labelcolor=seq_time_color)

    ## Plot average line
    ax2.hlines([1], -1, len(individual_results) + 1, linewidth=0.8, linestyles="dotted")
    # ax2.hlines([mean], -1, len(individual_results) + 1, linewidth=1.2)
    mean_xs = range(-1, len(individual_results) + 1)
    ax2.plot(mean_xs, [mean for _ in mean_xs], linewidth=1.2, label='Mean Speedup')

    ax1.yaxis.tick_right()
    ax1.yaxis.set_label_position("right")
    ax2.yaxis.tick_left()
    ax2.yaxis.set_label_position("left")
    old_ticks = ax2.get_yticks()
    # print(old_ticks)
    ax2.set_yticks([1] + old_ticks[:-1])
    # plt.yscale("log")
    # plt.yticks(range(1, 18, 2))
    # plt.ylim((0.1, 20))
    ax2.legend()
    # ax2.grid()
    # plt.title("Unix50 Individual Speedups")
    plt.tight_layout()
    plt.savefig(os.path.join('../evaluation/plots', "unix50{}_individual_speedups_{}.pdf".format(small_prefix, parallelism)),bbox_inches='tight')

def make_coarse_unix50_bar_chart(all_results, scaleup_numbers, parallelism):
    ## Plot individual speedups
    individual_results = [distr_exec_speedup[scaleup_numbers.index(parallelism)]
                          for distr_exec_speedup, _ in all_results]
    absolute_seq_times_s = [absolute_seq / 1000
                            for _, absolute_seq in all_results]
    print("Unix50 individual speedups for {} parallelism:".format(parallelism), individual_results)
    
    aggrs = compute_and_print_aggrs(individual_results, absolute_seq_times_s)
    mean = aggrs[0]

    w = 0.2
    ind = np.arange(len(individual_results))
    speedup_color = 'tab:blue'
    seq_time_color = 'tab:red'

    fig, ax1 = plt.subplots()
    # print(fig.get_size_inches())
    fig.set_size_inches(10, 4)

    ## Plot speedup
    ax1.set_ylabel('Speedup', color=speedup_color)
    ax1.set_xlabel('Pipeline')
    plt.xlim(-1, len(individual_results))
    plt.hlines([1], -1, len(individual_results) + 1, linewidth=0.8)
    ax1.bar(ind-w, individual_results, width=2*w, align='center', color=speedup_color)
    ax1.tick_params(axis='y', labelcolor=speedup_color)

    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
    ax2.set_ylabel('Sequential Time (s)', color=seq_time_color)
    ax2.set_yscale("log")
    ax2.bar(ind+w, absolute_seq_times_s, width=2*w, align='center', color=seq_time_color)
    ax2.tick_params(axis='y', labelcolor=seq_time_color)
    # plt.yscale("log")
    # plt.yticks(range(1, 18, 2))
    # plt.ylim((0.1, 20))
    # plt.legend(loc='lower right')
    # plt.title("Unix50 Individual Speedups")
    plt.tight_layout()
    plt.savefig(os.path.join('../evaluation/plots', "unix50_coarse_individual_speedups_{}.pdf".format(parallelism)),bbox_inches='tight')

def get_pipelines_res(individual_results, absolute_seq_times_s, pipelines):
    speedups = [individual_results[index] for index in pipelines]
    seq_times = [absolute_seq_times_s[index] for index in pipelines]
    return speedups, seq_times

def make_unix50_scatter_plot(all_results, scaleup_numbers, parallelism):
    individual_results = [distr_exec_speedup[scaleup_numbers.index(parallelism)]
                          for distr_exec_speedup, _ in all_results]
    absolute_seq_times_s = [absolute_seq / 1000
                            for _, absolute_seq in all_results]
    print("Unix50 individual speedups for {} scripts and {} parallelism:".format(len(individual_results), parallelism), individual_results)
    
    aggrs = compute_and_print_aggrs(individual_results, absolute_seq_times_s)
    mean = aggrs[0]

    slowdown_pipelines = [2,19,31]
    no_speedup_pipelines = [13, 24, 25, 26, 29, 30]
    sort_pipelines = [0, 1, 3, 15, 16, 18, 20, 27, 28, 33]
    deep_pipelines = [7, 8, 9, 10, 11, 12, 14, 27, 28]
    io_pipelines = [4, 5, 6, 7, 8, 12, 14, 22, 23]
    special_pipelines = slowdown_pipelines + no_speedup_pipelines + sort_pipelines + deep_pipelines + io_pipelines
    rest_pipelines_set = set(range(len(all_results))) - set(special_pipelines)
    rest_pipelines = sorted(list(rest_pipelines_set))


    fig, ax = plt.subplots()
    fig.set_size_inches(10, 5)

    ## Slowdown Pipelines
    speeds, seqs = get_pipelines_res(individual_results, absolute_seq_times_s, rest_pipelines)
    ax.scatter(seqs, speeds, label="Highly Parallelizable")
    print("|-- Highly parallelizable speedups:", speeds)

    speeds, seqs = get_pipelines_res(individual_results, absolute_seq_times_s, sort_pipelines)
    ax.scatter(seqs, speeds, label="Contain sort")

    speeds, seqs = get_pipelines_res(individual_results, absolute_seq_times_s, io_pipelines)
    ax.scatter(seqs, speeds, label="Non CPU-heavy")

    speeds, seqs = get_pipelines_res(individual_results, absolute_seq_times_s, deep_pipelines)
    ax.scatter(seqs, speeds, label="Deep")

    speeds, seqs = get_pipelines_res(individual_results, absolute_seq_times_s, no_speedup_pipelines)
    ax.scatter(seqs, speeds, label="Non parallelizable")
    print("|-- Non parallelizable speedups:", speeds)

    speeds, seqs = get_pipelines_res(individual_results, absolute_seq_times_s, slowdown_pipelines)
    ax.scatter(seqs, speeds, label="Contain head")
    print("|-- Head script speedups:", speeds)

    ax.grid()
    ax.set_ylabel('Speedup')
    ax.set_xlabel('Sequential Time (s)')
    ax.set_xscale("log")
    # plt.hlines([1], -1, 1000000, linewidth=0.8)
    old_ylim = plt.ylim()
    plt.yticks(list(plt.yticks()[0]) + [1])
    plt.xlim(0, max(absolute_seq_times_s) * 2)
    plt.ylim(old_ylim)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join('../evaluation/plots', "unix50_scatter_speedups_{}.pdf".format(parallelism)),bbox_inches='tight')

def plot_unix50_avg_speedup(all_results, scaleup_numbers, filename):
    avg_results = aggregate_unix50_results(all_results, scaleup_numbers)
    print("Unix50 average speedup:", avg_results)

    ## Plot average speedup
    fig, ax = plt.subplots()

    ## Plot speedup
    ax.set_ylabel('Speedup')
    ax.set_xlabel('--width')
    ax.plot(scaleup_numbers, avg_results, '-o', linewidth=0.5, label='Parallel')
    # plt.yscale("log")
    plt.xticks(scaleup_numbers)
    plt.legend(loc='lower right')
    plt.title("Unix50 Throughput")
    plt.tight_layout()
    plt.savefig(os.path.join('../evaluation/plots', filename),bbox_inches='tight')

def collect_all_unix50_results(unix50_results_dir, scaleup_numbers=[2, 4, 8, 16], suffix='distr.time'):

    files = [f for f in os.listdir(unix50_results_dir)]
    # print(files)
    pipeline_numbers = sorted(list(set([f.split('_')[2] for f in files])))
    # print(pipeline_numbers)

    all_results = [collect_unix50_pipeline_scaleup_times(pipeline_number,
                                                         unix50_results_dir,
                                                         scaleup_numbers,
                                                         suffix=suffix)
                   for pipeline_number in pipeline_numbers]
    fan_in_fan_out_results = [collect_unix50_pipeline_fan_in_fan_out_scaleup_times(pipeline_number,
                                                                    unix50_results_dir,
                                                                    scaleup_numbers)
                   for pipeline_number in pipeline_numbers]
    
    return (all_results, fan_in_fan_out_results)

def collect_unix50_scaleup_times(all_results, scaleup=[2,4,8,16], small_prefix=""):
    
    # print(all_results)

    for parallelism in scaleup:
        make_unix50_bar_chart(all_results, scaleup, parallelism, small_prefix=small_prefix)
        make_unix50_scatter_plot(all_results, scaleup, parallelism)

    plot_unix50_avg_speedup(all_results, scaleup, "unix50_throughput_scaleup.pdf")
    

# def collect_unix50_coarse_scaleup_times(all_results):
#     print(len(all_results))
#     scaleup_numbers = [16]

#     # print(all_results)

#     for parallelism in scaleup_numbers:
#         make_coarse_unix50_bar_chart(all_results, scaleup_numbers, parallelism)

#     plot_unix50_avg_speedup(all_results, scaleup_numbers, "unix50_coarse_throughput_scaleup.pdf")


def get_statistics_from_lines(lines):
    minimums = lines[0].get_ydata().copy()
    maximums = lines[0].get_ydata().copy()
    sums = [0 for _ in lines[0].get_ydata()]
    for line in lines:
        for i, y in enumerate(line.get_ydata()):
            minimums[i] = min(minimums[i], y)
            maximums[i] = max(maximums[i], y)
            sums[i] += y
    avgs = [float(s)/len(lines) for s in sums]
    return minimums, maximums, avgs

def set_tiling_axes_labels_ticks(fig):
    axs = fig.get_axes()
    for ax in axs:
        if(ax.is_first_col()):
            ax.set_ylabel('Speedup')
        if(ax.is_last_row()):
            ax.set_xlabel('--width')
        if(not ax.is_last_row()):
            ax.set_xticklabels([])
        # ax.label_outer()

def plot_tiling_experiments(fig, gs, experiments, all_experiment_results, all_scaleup_numbers,
                            correctness=None, custom_scaleup_plots={},
                            line_plot_configs=default_line_plot_configs):
    total_lines = []
    averages = [[] for _ in all_scaleup_numbers]
    no_eager_averages = [[] for _ in all_scaleup_numbers]
    ## Plot microbenchmarks
    for i, experiment in enumerate(experiments):
        ax = fig.add_subplot(gs[i])
        all_speedup_results = all_experiment_results[experiment]
        lines, best_result, no_eager_result = plot_scaleup_lines(experiment, all_scaleup_numbers, all_speedup_results, 
                                                                 custom_scaleup_plots, ax, line_plot_configs=line_plot_configs)
        # if(experiment == "double_sort"):
        total_lines += lines
        ax.set_xticks(all_scaleup_numbers[1:])

        text_color = 'black'
        if (not correctness is None 
            and any_wrong(correctness, experiment, all_line_plots)):
            text_color = 'red'
        ax.text(.5,.91,pretty_names[experiment],
                horizontalalignment='center',
                transform=ax.transAxes,
                color=text_color)
        # ax.set_yticks([])
        fig.add_subplot(ax)

        ## Update averages
        for i, res in enumerate(best_result):
            averages[i].append(res)
        for i, res in enumerate(no_eager_result):
            no_eager_averages[i].append(res)
    
    set_tiling_axes_labels_ticks(fig)
    
    return (total_lines, averages, no_eager_averages)


def print_aggregates(prefix, averages, no_eager_averages):
    ## Print average, geo-mean
    one_liner_averages = [sum(res)/len(res) for res in averages]
    all_no_eager_averages = [sum(res)/len(res) for res in no_eager_averages]
    geo_means = [math.exp(np.log(res).sum() / len(res))
                 for res in averages]
    print(prefix, "One-liners Aggregated results:")
    print(" |-- Averages:", one_liner_averages)
    print(" |-- No Eager Averages:", all_no_eager_averages)
    print(" |-- Geometric Means:", geo_means)

## TODO: Add more experiments to be ploted in the report
def report_all_one_liners(all_scaleup_numbers, all_experiment_results, correctness):
    global all_line_plots

    line_plots = ["split", "eager", "blocking-eager", "no-eager", "no-aux-cat-split"]
    legend_names = ["PaSh",
                    "PaSh w/o split",
                    "Blocking Eager",
                    "No Eager",
                    "No Aux Cat-Split"]

    fig = plt.figure()
    columns = 5
    rows = (len(all_experiment_results.keys()) // columns) + 1
    gs = fig.add_gridspec(rows, columns, hspace=0.05)

    plot_res = plot_tiling_experiments(fig, gs, all_experiments, 
                                       all_experiment_results, 
                                       all_scaleup_numbers, 
                                       correctness)
    _total_lines, averages, no_eager_averages = plot_res 

    legend_artists = [default_line_plot_configs[lplot].get_artist() for lplot in line_plots]
    plt.legend(legend_artists, legend_names, loc='lower right', fontsize=16)
    # plt.title(pretty_names[experiment])

    fig.set_size_inches(columns * 6, rows * 5)
    plt.tight_layout()
    plt.savefig(os.path.join('../evaluation/plots', "all_one_liners_scaleup.pdf"),bbox_inches='tight')

    print_aggregates("All", averages, no_eager_averages)

def plot_one_liners_tiling(all_experiment_results, experiments,
                           custom_scaleup_plots,
                           all_scaleup_numbers=[2, 4, 8, 16, 32, 64],
                           prefix=""):

    line_plots = ["split", "eager", "blocking-eager", "no-eager"]

    legend_names = ["PaSh",
                    "PaSh w/o split",
                    "Blocking Eager",
                    "No Eager"]

    fig = plt.figure()
    gs = fig.add_gridspec(2, 5, hspace=0.05)
    
    plot_res = plot_tiling_experiments(fig, gs, experiments, 
                                       all_experiment_results, 
                                       all_scaleup_numbers,
                                       custom_scaleup_plots=custom_scaleup_plots)
    _total_lines, averages, no_eager_averages = plot_res

    legend_artists = [default_line_plot_configs[lplot].get_artist() for lplot in line_plots]
    plt.legend(legend_artists, legend_names, loc='lower right', fontsize=16)
    # plt.title(pretty_names[experiment])

    fig.set_size_inches(27, 8.2)
    plt.tight_layout()
    plt.savefig(os.path.join('../evaluation/plots', "{}tiling_throughput_scaleup.pdf".format(prefix)),bbox_inches='tight')

    print_aggregates("Systems", averages, no_eager_averages)

def plot_less_one_liners_tiling(all_experiment_results, all_sequential_results, experiments, unix50_results, oasa_results):

    all_scaleup_numbers = [2, 4, 8, 16, 32, 64]

    coarse_custom_scaleup_plots = {"minimal_grep" : ["blocking-eager"],
                                   "minimal_sort" : ["eager"],
                                   "topn" : ["no-eager"],
                                   "wf" : ["no-eager"],
                                   "spell" : ["mini-split"],
                                   "diff" : ["no-eager"],
                                   "bigrams" : ["mini-split"],
                                   "set-diff" : ["eager"],
                                   "shortest_scripts" : ["no-eager"],
                                  }

    line_plot_configs = {'eager': LinePlotConfig('-o', 'tab:blue', 'Parallel w/o split'),
                         'split': LinePlotConfig('-o', 'tab:blue', 'Parallel'),
                         'mini-split': LinePlotConfig('-o', 'tab:blue', 'Parallel'),
                         'blocking-eager': LinePlotConfig('-o', 'tab:blue', 'Blocking Eager'),
                         'no-eager': LinePlotConfig('-o', 'tab:blue', 'No Eager')}

    # confs = ["Parallel"]

    fig = plt.figure()
    gs = fig.add_gridspec(3, 3, hspace=0.05)

    plot_res = plot_tiling_experiments(fig, gs, experiments, 
                                       all_experiment_results, 
                                       all_scaleup_numbers,
                                       custom_scaleup_plots=coarse_custom_scaleup_plots,
                                       line_plot_configs=line_plot_configs)
    total_lines, averages, no_eager_averages = plot_res

    # plt.legend(total_lines, confs, loc='lower right', fontsize=16)
    # plt.title(pretty_names[experiment])

    fig.set_size_inches(22.5, 12.5)
    plt.tight_layout()
    ## TODO: Replace the prefix with a constant
    plt.savefig(os.path.join('../evaluation/plots', "coarse_tiling_throughput_scaleup.pdf"),bbox_inches='tight')

    ## Plot the aggregate plot for the one liners
    plot_less_one_liners_aggregate(total_lines, all_scaleup_numbers)
    plot_all_less_one_liners_one_plot(total_lines, all_scaleup_numbers, experiments)

    ## Plot two bars, one for the fan-in fan-out and one for the total
    scaleup_number = 16
    plot_bar_chart_one_liners(total_lines, all_experiment_results, all_sequential_results, experiments, scaleup_number)
    plot_bar_chart_one_liners_and_unix50(total_lines, all_experiment_results, all_sequential_results, experiments, scaleup_number, unix50_results, oasa_results)
    plot_one_bar_chart_one_liners_and_unix50(total_lines, all_experiment_results, all_sequential_results, experiments, scaleup_number, unix50_results)

    print_aggregates("Coarse", averages, no_eager_averages)

def gather_abs_times_speedups_bar(total_lines, all_experiment_results, all_sequential_results, experiments, scaleup_number):
    seq_results_s = [all_sequential_results[experiment] / 1000.0 for experiment in experiments]

    ## Gather the good results
    good_speedups = []
    for l_i, line in enumerate(total_lines):
        i, = np.where(line.get_xdata() == scaleup_number)
        speedup = line.get_ydata()[int(i)]
        good_speedups.append(speedup)
    # print(good_results)
    good_results = [seq_results_s[l_i] / speedup for i, speedup in enumerate(good_speedups)]

    ## Gather the no-aux transformation results
    no_aux_speedups = []
    for ex_i, experiment in enumerate(experiments):
        line = total_lines[ex_i]
        no_aux_res = all_experiment_results[experiment]["no-aux-cat-split"]
        i, = np.where(line.get_xdata() == scaleup_number) 
        speedup = no_aux_res[int(i)]
        no_aux_speedups.append(speedup)
    # print(no_aux_results)

    no_aux_results = [safe_zero_div(seq_results_s[i], speedup)
                      for i, speedup in enumerate(no_aux_speedups)]
    return (good_results, good_speedups, no_aux_results, no_aux_speedups, seq_results_s)


def plot_bar_chart_one_liners(total_lines, all_experiment_results, all_sequential_results, experiments, scaleup_number):
    res = gather_abs_times_speedups_bar(total_lines, all_experiment_results, all_sequential_results, experiments, scaleup_number)
    (good_results, good_speedups, no_aux_results, no_aux_speedups, seq_results_s) = res
    w = 0.12
    ind = np.arange(len(good_results))
    good_speedup_color = 'tab:blue'
    no_aux_speedup_color = 'tab:orange'
    seq_color = 'tab:green'

    fig, ax = plt.subplots()
    # print(fig.get_size_inches())
    fig.set_size_inches(12, 6)

    ## Plot speedup
    ax.set_xlabel('Execution time (s)')
    ax.set_ylabel('Script')
    ax.grid(axis='x', zorder=0)

    # plt.vlines([1], -1, len(good_results) + 1, linewidth=0.8)
    ax.barh(ind+2*w, good_results[::-1], height=2*w, align='center', 
            color=good_speedup_color, label='Par', zorder=3)
    ax.barh(ind, no_aux_results[::-1], height=2*w, align='center', 
            color=no_aux_speedup_color, label='Par -aux', zorder=3)
    ax.barh(ind-2*w, seq_results_s[::-1], height=2*w, align='center', 
            color=seq_color, label='Seq', zorder=3)
    ylabels = [pretty_names[exp] for exp in experiments]
    plt.yticks(ind, ylabels[::-1])    
    # plt.ylim((0.1, 20))
    plt.legend(loc='lower right')
    plt.tight_layout()

    plt.savefig(os.path.join('../evaluation/plots', "coarse_one_liners_bar_{}.pdf".format(scaleup_number)),bbox_inches='tight')
    print("Transformation averages:")
    print("|-- All transformations:", sum(good_speedups) / len(good_speedups))
    print("|-- No Cat-Split transformation:", sum(no_aux_speedups) / len(no_aux_speedups))


def plot_bar_chart_one_liners_and_unix50(total_lines, all_experiment_results, all_sequential_results, experiments, scaleup_number, unix50_all_results, oasa_results):
    res = gather_abs_times_speedups_bar(total_lines, all_experiment_results, all_sequential_results, experiments, scaleup_number)
    (good_results, good_speedups, no_aux_results, no_aux_speedups, seq_results_s) = res

    unix50_results, unix50_results_fan_in = unix50_all_results
    index_16 = 3
    unix50_seq_times_s = [absolute_seq / 1000 for _, absolute_seq in unix50_results]
    unix50_good_speedups = [float(distr_exec_speedup[index_16])
                            for distr_exec_speedup, _ in unix50_results]
    unix50_good_results = [unix50_seq_times_s[i] / float(distr_exec_speedup)
                           for i, distr_exec_speedup in enumerate(unix50_good_speedups)]
    unix50_no_aux_speedups = [distr_exec_speedup[index_16]
                              for distr_exec_speedup, _ in unix50_results_fan_in]
    unix50_no_aux_results = [safe_zero_div(unix50_seq_times_s[i], distr_exec_speedup)
                             for i, distr_exec_speedup in enumerate(unix50_no_aux_speedups)]

    ## Extract the OASA experiment data
    ##
    ## TODO: Disgusting copy-paste...
    oasa_experiments = [str(i) for i in range(1,5)]
    oasa_seq_times_s = [oasa_results[exp]['seq'] / 1000.0 
                       for exp in oasa_experiments]
    oasa_good_speedups = [float(oasa_results[exp]['split'][index_16]) 
                          for exp in oasa_experiments]
    oasa_good_results = [safe_zero_div(oasa_seq_times_s[i], float(distr_exec_speedup))
                         for i, distr_exec_speedup in enumerate(oasa_good_speedups)]
    oasa_no_aux_speedups = [float(oasa_results[exp]['no-aux-cat-split'][index_16]) 
                            for exp in oasa_experiments]
    oasa_no_aux_results = [safe_zero_div(oasa_seq_times_s[i], distr_exec_speedup)
                             for i, distr_exec_speedup in enumerate(oasa_no_aux_speedups)]

    all_good_results = good_results + unix50_good_results + oasa_good_results
    all_no_aux_results = no_aux_results + unix50_no_aux_results + oasa_no_aux_results
    all_seq_results = seq_results_s + unix50_seq_times_s + oasa_seq_times_s

    ylim = 1500
    w = 0.12
    ind = np.arange(len(all_good_results))
    good_speedup_color = 'tab:blue'
    no_aux_speedup_color = 'tab:orange'
    seq_color = 'tab:green'

    fig, ax = plt.subplots()
    # print(fig.get_size_inches())
    fig.set_size_inches(25, 7.5)

    ## Plot speedup
    ax.set_ylabel('Execution time (s)')
    ax.set_xlabel('Script')
    ax.grid(axis='y', zorder=0)

    # plt.vlines([1], -1, len(good_results) + 1, linewidth=0.8)
    ax.bar(ind+2*w, all_good_results, width=2*w, align='center', 
            color=good_speedup_color, label='Parallel', zorder=3)
    ax.bar(ind, all_no_aux_results, width=2*w, align='center', 
            color=no_aux_speedup_color, label='No Cat-Split', zorder=3)
    ax.bar(ind-2*w, all_seq_results, width=2*w, align='center', 
            color=seq_color, label='Baseline', zorder=3)
    ## Add text on top
    prev_i = -2
    for i, v in enumerate(all_seq_results):
        if(v > ylim):
            if(prev_i == i-1):
                y_pos = ylim*1.07
            else:
                y_pos = ylim * 1.01
                prev_i = i
            plt.text(i-2*w, y_pos, str(int(v)), ha='center')

    xlabels = [pretty_names[exp] for exp in experiments] + ["unix50-{}".format(i) for i in range(len(unix50_good_results))] + ["buses-{}".format(i) for i in range(1,5)]
    plt.xticks(ind, xlabels, rotation=45, ha="right") 
    # ax.set_xscale("log")
    plt.ylim((0, ylim))
    plt.legend(loc='upper right')
    plt.tight_layout()

    ## Add vertical lines
    line_x=len(good_results) - 4 * w
    plt.axvline(x=line_x, color='black', linestyle='--')

    line_x=len(good_results + unix50_good_results) - 4 * w
    plt.axvline(x=line_x, color='black', linestyle='--')

    all_good_speedups = good_speedups + unix50_good_speedups + oasa_good_speedups
    all_no_aux_speedups = no_aux_speedups + unix50_no_aux_speedups + oasa_no_aux_speedups
    print(all_no_aux_speedups)
    print("All Transformation averages:")
    print("|-- All transformations:", sum(all_good_speedups) / len(all_good_speedups))
    print("|-- No Cat-Split transformation:", sum(all_no_aux_speedups) / len(all_no_aux_speedups))

    print("All Transformation minimum:")
    print("|-- All transformations:", min(all_good_speedups))
    print("|-- No Cat-Split transformation:", min(all_no_aux_speedups))

    print("All Transformation maximum:")
    print("|-- All transformations:", max(all_good_speedups))
    print("|-- No Cat-Split transformation:", max(all_no_aux_speedups))


    plt.savefig(os.path.join('../evaluation/plots', "coarse_all_bar_{}.pdf".format(scaleup_number)),bbox_inches='tight')

def plot_one_bar_chart_one_liners_and_unix50(total_lines, all_experiment_results, all_sequential_results, experiments, scaleup_number, unix50_all_results):
    res = gather_abs_times_speedups_bar(total_lines, all_experiment_results, all_sequential_results, experiments, scaleup_number)
    (_, good_speedups, _, _, seq_results_s) = res

    # nfa_index = experiments.index('minimal_grep')
    # good_speedups[nfa_index] = 15.81

    unix50_results, _ = unix50_all_results
    index_16 = 3
    unix50_good_speedups = [float(distr_exec_speedup[index_16])
                            for distr_exec_speedup, _ in unix50_results]
    unix50_good_speedups = [res for res in unix50_good_speedups
                            if res > 0.2]

    w = 0.2
    ind = np.arange(len(good_speedups + unix50_good_speedups))
    good_speedup_color = 'tab:blue'

    fig, ax = plt.subplots()
    # print(fig.get_size_inches())
    fig.set_size_inches(25, 7.5)

    ## Plot speedup
    ax.set_ylabel('Speedup')
    ax.set_xlabel('Script')
    ax.grid(axis='y', zorder=0)

    ax.bar(ind, good_speedups + unix50_good_speedups, width=2*w, align='center', 
            color=good_speedup_color, label='PaSh', zorder=3)

    xlabels = [pretty_names[exp] for exp in experiments] + ["unix50-{}".format(i) for i in range(len(unix50_good_speedups))]
    plt.xticks(ind, xlabels, rotation=45, ha="right")
    old_ylim = plt.ylim()
    plt.yticks(list(plt.yticks()[0]) + [1])
    plt.ylim(old_ylim)
    # ax.set_xscale("log")
    # plt.legend(loc='upper right')
    plt.tight_layout()

    all_good_speedups = good_speedups + unix50_good_speedups
    print("One bar all Transformation averages:")
    print("|-- All transformations:", sum(all_good_speedups) / len(all_good_speedups))
    print("One bar all Transformation minimum:")
    print("|-- All transformations:", min(all_good_speedups))
    print("One bar all Transformation maximum:")
    print("|-- All transformations:", max(all_good_speedups))


    plt.savefig(os.path.join('../evaluation/plots', "coarse_all_one_bar_{}.pdf".format(scaleup_number)),bbox_inches='tight')



def plot_less_one_liners_aggregate(total_lines, all_scaleup_numbers):
    mins, maxs, avgs = get_statistics_from_lines(total_lines)

    ## Plot one plot with all together
    fig, ax = plt.subplots()
    ax.plot(all_scaleup_numbers, avgs)
    ax.fill_between(all_scaleup_numbers, mins, maxs, alpha=0.2)
    ax.hlines([1], all_scaleup_numbers[0], all_scaleup_numbers[-1], linewidth=0.8, linestyles="dotted")
    loc = plticker.MultipleLocator(base=10.0) # this locator puts ticks at regular intervals
    ax.yaxis.set_major_locator(loc)
    plt.xticks(all_scaleup_numbers)
    plt.xlim(all_scaleup_numbers[0], all_scaleup_numbers[-1])
    ax.grid()
    plt.ylim(0,maxs[-1]*1.05)
    fig.set_size_inches(9, 5)
    set_tiling_axes_labels_ticks(fig)
    plt.tight_layout()
    plt.savefig(os.path.join('../evaluation/plots', "coarse_aggregate_throughput_scaleup.pdf"),bbox_inches='tight')

def plot_all_less_one_liners_one_plot(total_lines, all_scaleup_numbers, experiments):
    mins, maxs, avgs = get_statistics_from_lines(total_lines)
    fig, ax = plt.subplots()
    ax.plot(all_scaleup_numbers, avgs, label="Average", linewidth=2)
    for i, line in enumerate(total_lines):
        label = pretty_names[experiments[i]]
        ax.plot(line.get_xdata(), line.get_ydata(), alpha=0.8, label=label, linewidth=1, linestyle="dashed")
        # new_line = pltlines.Line2D(line.get_xdata(), line.get_ydata())
        # new_line.update_from(line)
        # ax.add_line(new_line)
    ax.hlines([1], all_scaleup_numbers[0], all_scaleup_numbers[-1], linewidth=0.8, linestyles="dotted")
    loc = plticker.MultipleLocator(base=10.0) # this locator puts ticks at regular intervals
    ax.yaxis.set_major_locator(loc)
    plt.xticks(all_scaleup_numbers)
    plt.xlim(all_scaleup_numbers[0], all_scaleup_numbers[-1])
    ax.grid()
    # ax.set_yscale("log")
    plt.ylim(0,maxs[-1]*1.05)
    # plt.ylim(0,20)
    fig.set_size_inches(9, 5)
    set_tiling_axes_labels_ticks(fig)
    plt.tight_layout()
    plt.legend(loc='upper left', fontsize=14)
    plt.savefig(os.path.join('../evaluation/plots', "coarse_all_in_one_throughput_scaleup.pdf"),bbox_inches='tight')


def format_correctness(correctness):
    global all_experiments
    global all_line_plots

    print("\nWrong Result Summary (If results are missing for complete conf they are ignored):")
    for experiment in all_experiments:
        for line_plot in all_line_plots:
            result_correctness = correctness[experiment][line_plot]
            all_missing = all([res == "Not Found" for n, res in result_correctness])
            wrong_ns = []
            missing_ns = []
            for n, res in result_correctness:
                if res == "Wrong":
                    wrong_ns.append(n)
                elif (res == "Not Found"):
                    missing_ns.append(n)
            if len(wrong_ns) > 0:
                print("|-- WARNING: Wrong output for", experiment, "-", line_plot, "-", wrong_ns)
            if len(missing_ns) > 0 and not all_missing:
                print("|-- WARNING: Missing output for", experiment, "-", line_plot, "-", missing_ns)            

def any_wrong(correctness, experiment, line_plots):
    for line_plot in line_plots:
        result_correctness = correctness[experiment][line_plot]
        any_wrong = any([res == "Wrong" for n, res in result_correctness])
        if any_wrong:
            return True
    return False


## Set the fonts to be larger
SMALL_SIZE = 22
MEDIUM_SIZE = 24
BIGGER_SIZE = 30

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

# Plot microbenchmarks
diff_results = []
all_scaleup_numbers = [2, 4, 8, 16, 32, 64]
all_experiment_results = {}
all_sequential_results = {}
correctness = {}
for experiment in all_experiments:
    all_speedup_results, output_diff, sequential_time = collect_scaleup_line_speedups(experiment, all_scaleup_numbers, RESULTS)
    all_experiment_results[experiment] = all_speedup_results
    all_sequential_results[experiment] = sequential_time
    correctness[experiment] = output_diff

small_one_liners_scaleup_numbers = [2, 16]
small_one_liner_results = {}
for experiment in all_experiments:
    small_speedup_results, _, sequential_time = collect_scaleup_line_speedups(experiment, small_one_liners_scaleup_numbers, SMALL_RESULTS)
    small_one_liner_results[experiment] = small_speedup_results

## Gather results for OASA benchmarks
oasa_experiments = [str(i) for i in range(1,5)]
oasa_scaleup_numbers = all_scaleup_numbers
oasa_experiment_results = {}
for experiment in oasa_experiments:
    all_speedup_results, _, sequential_time = collect_scaleup_line_speedups(experiment, oasa_scaleup_numbers, OASA_RESULTS)
    oasa_experiment_results[experiment] = all_speedup_results
    oasa_experiment_results[experiment]['seq'] = sequential_time

## Make a report of all one-liners
report_all_one_liners(all_scaleup_numbers, all_experiment_results, correctness)

## Legacy unix50 results
if args.all:
    unix50_results, unix50_results_fan_in = collect_all_unix50_results(UNIX50_RESULTS)

## Collect all unix50 results
if args.eurosys2021:
    small_unix50_results, _ = collect_all_unix50_results(SMALL_UNIX50_RESULTS, scaleup_numbers=[4], suffix='distr_auto_split.time')
    big_unix50_results, _ = collect_all_unix50_results(BIG_UNIX50_RESULTS, scaleup_numbers=[16], suffix='distr_auto_split.time')

##
## Theory Paper
##
coarse_experiments = ["minimal_grep",
                      "minimal_sort",
                      "topn",
                      "wf",
                      "spell",
                      "bigrams",
                      "diff",
                      "set-diff",
                      "shortest_scripts"]

if args.all:
    plot_less_one_liners_tiling(all_experiment_results, all_sequential_results, 
                                coarse_experiments, (unix50_results, unix50_results_fan_in),
                                oasa_experiment_results)
    generate_tex_coarse_table(coarse_experiments)
    # collect_unix50_coarse_scaleup_times(unix50_results)


##
## Systems Paper
##
experiments = ["minimal_grep",
               "minimal_sort",
               "topn",
               "wf",
               "spell",
               "diff",
               "bigrams",
               "set-diff",
               "double_sort",
               "shortest_scripts"]

## Large inputs `-l`
custom_scaleup_plots = {"minimal_grep" : ["eager", "blocking-eager"],
                        "minimal_sort": ["eager", "blocking-eager", "no-eager"],
                        "topn": ["eager", "blocking-eager", "no-eager"],
                        "wf": ["eager", "blocking-eager", "no-eager"],
                        "spell" : ["split", "eager"],
                        "diff" : ["eager", "blocking-eager", "no-eager"],
                        "bigrams" : ["split", "eager"],
                        "set-diff" : ["eager", "blocking-eager", "no-eager"],
                        "double_sort" : ["split", "eager", "blocking-eager", "no-eager"],
                        "shortest_scripts" : ["eager", "blocking-eager", "no-eager"]}

if args.eurosys2021:
    plot_one_liners_tiling(all_experiment_results, experiments, custom_scaleup_plots)

## Medium input `-m`
medium_custom_scaleup_plots = {"minimal_grep" : ["split", "blocking-eager"],
                               "minimal_sort": ["split", "blocking-eager", "no-eager"],
                               "topn": ["split", "blocking-eager", "no-eager"],
                               "wf": ["split", "blocking-eager", "no-eager"],
                               "spell" : ["split", "eager"],
                               "diff" : ["split", "blocking-eager", "no-eager"],
                               "bigrams" : ["split", "eager"],
                               "set-diff" : ["split", "blocking-eager", "no-eager"],
                               "double_sort" : ["split", "eager", "blocking-eager", "no-eager"],
                               "shortest_scripts" : ["split", "blocking-eager", "no-eager"]}


if args.eurosys2021:
    plot_one_liners_tiling(small_one_liner_results, experiments,
                           medium_custom_scaleup_plots,
                           all_scaleup_numbers=small_one_liners_scaleup_numbers,
                           prefix="medium_")

## Small input `-s`
small_custom_scaleup_plots = {"minimal_grep" : ["split"],
                              "minimal_sort": ["split"],
                              "topn": ["split"],
                              "wf": ["split"],
                              "spell" : ["split"],
                              "diff" : ["split"],
                              "bigrams" : ["split"],
                              "set-diff" : ["split"],
                              "double_sort" : ["split"],
                              "shortest_scripts" : ["split"]}

if args.eurosys2021:
    plot_one_liners_tiling(small_one_liner_results, experiments,
                           small_custom_scaleup_plots,
                           all_scaleup_numbers=small_one_liners_scaleup_numbers,
                           prefix="small_")

if args.eurosys2021:
    generate_tables(experiments, results_dir=SMALL_RESULTS, table_suffix="-small", small=True)
    generate_tables(experiments)
    collect_unix50_scaleup_times(small_unix50_results, scaleup=[4], small_prefix="_1GB")
    collect_unix50_scaleup_times(big_unix50_results, scaleup=[16], small_prefix="_10GB")
    plot_sort_with_baseline(RESULTS)

## Legacy plots
if args.all:
    collect_unix50_scaleup_times(unix50_results)


## Format and print correctness results
if args.all:
    format_correctness(correctness)
