import sys
import math
import os
import numpy as np
import matplotlib.pyplot as plt
import statistics

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
UNIX50_RESULTS = "../evaluation/results/unix50/"
COARSE_UNIX50_RESULTS = "../evaluation/results/unix50-naive/"


experiments = ["minimal_grep",
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

coarse_experiments = ["minimal_grep",
                      "minimal_sort",
                      "topn",
                      "wf",
                      "spell",
                      "diff",
                      "bigrams",
                      "set-diff",
                      "shortest_scripts"]


pretty_names = {"minimal_grep" : "grep",
                "minimal_sort" : "sort",
                "wf" : "wf",
                "topn" : "top-n",
                "grep" : "grep-light",
                "bigrams" : "bi-grams",
                "alt_bigrams" : "bi-grams-opt",
                "spell" : "spell",
                "shortest_scripts" : "shortest-scripts",
                "diff" : "diff",
                "set-diff" : "set-diff",
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

custom_scaleup_plots = {"set-diff" : ["eager", "blocking-eager", "no-eager"],
                        "spell" : ["split", "eager"],
                        "bigrams" : ["split", "eager"]}

coarse_custom_scaleup_plots = {"minimal_grep" : ["blocking-eager"],
                               "minimal_sort" : ["no-eager"],
                               "topn" : ["no-eager"],
                               "wf" : ["no-eager"],
                               "spell" : ["mini-split"],
                               "diff" : ["no-eager"],
                               "bigrams" : ["mini-split"],
                               "set-diff" : ["no-eager"],
                               "shortest_scripts" : ["no-eager"],
                               }

input_filename_sizes = {"1G": "1~GB",
                        "3G": "3~GB",
                        "10G": "10~GB",
                        "100G": "100~GB",
                        "1M": "1~MB",
                        "10M": "10~MB",
                        "100M": "100~MB",
                        "all_cmds_x1000": "85~MB"}

def safe_zero_div(a, b):
    if(b == 0):
        print("WARNING: Division by zero")
        return 0
    else:
        return a / b

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
            raise ValueError
        return sum(times)
    except:
        print("!! WARNING: Filename:", filename, "not found!!!")
        raise ValueError

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
                    return True
            return False
    except:
        print("!! WARNING: Filename:", filename, "not found!!!")
        return False

def collect_distr_experiment_execution_times(prefix, suffix, scaleup_numbers):
    numbers = [read_distr_execution_time('{}{}_{}'.format(prefix, n, suffix))
               for n in scaleup_numbers]
    return numbers

def collect_experiment_scaleup_times(prefix, scaleup_numbers):
    ## Since we have the same input size in all cases, only use the
    ## one sequential execution for the sequential time
    seq_numbers = [read_total_time('{}{}_seq.time'.format(prefix, scaleup_numbers[0]))
                   for _ in scaleup_numbers]
    distr_numbers = collect_distr_experiment_execution_times(prefix, 'distr.time', scaleup_numbers)
    compile_numbers = [read_distr_total_compilation_time('{}{}_distr.time'.format(prefix, n))
                       for n in scaleup_numbers]
    return (seq_numbers, distr_numbers, compile_numbers)

def collect_experiment_speedups(prefix, scaleup_numbers):
    seq_numbers, distr_numbers, compile_numbers = collect_experiment_scaleup_times(prefix, scaleup_numbers)
    distr_speedup = [safe_zero_div(seq_numbers[i], t) for i, t in enumerate(distr_numbers)]
    compile_distr_speedup = [safe_zero_div(seq_numbers[i], t + compile_numbers[i]) for i, t in enumerate(distr_numbers)]
    return (distr_speedup, compile_distr_speedup)

def collect_baseline_experiment_speedups(prefix, scaleup_numbers, base_seq):
    seq_numbers = [read_total_time('{}{}_seq.time'.format(prefix, n))
                   for n in scaleup_numbers]
    speedup = [safe_zero_div(base_seq, t) for t in seq_numbers]
    return speedup

def collect_distr_experiment_speedup(prefix, suffix, scaleup_numbers):
    # TODO: If this messed things up change it back
    # seq_numbers, _, _ = collect_experiment_scaleup_times(prefix, scaleup_numbers)
    seq_numbers = [read_total_time('{}{}_seq.time'.format(prefix, scaleup_numbers[0]))
                   for _ in scaleup_numbers]
    distr_numbers = collect_distr_experiment_execution_times(prefix, suffix, scaleup_numbers)
    distr_speedup = [safe_zero_div(seq_numbers[i], t) for i, t in enumerate(distr_numbers)]
    return distr_speedup

def collect_experiment_command_number(prefix, suffix, scaleup_numbers):
    command_numbers = [read_distr_command_number('{}{}_{}'.format(prefix, n, suffix))
                       for n in scaleup_numbers]
    return command_numbers

def check_output_diff_correctness(prefix, scaleup_numbers):
    wrong_diffs = [n for n in scaleup_numbers
                   if not check_output_diff_correctness_for_experiment('{}{}_distr.time'.format(prefix, n))]
    wrong_diffs_no_eager = [n for n in scaleup_numbers
                   if not check_output_diff_correctness_for_experiment('{}{}_distr_no_eager.time'.format(prefix, n))]
    wrong_diffs_no_task_par_eager = [n for n in scaleup_numbers
                   if not check_output_diff_correctness_for_experiment('{}{}_distr_no_task_par_eager.time'.format(prefix, n))]
    return (wrong_diffs, wrong_diffs_no_eager, wrong_diffs_no_task_par_eager)


def collect_scaleup_times_common(experiment, all_scaleup_numbers, results_dir, custom_scaleup_plots, ax):
    print(experiment)

    # all_scaleup_numbers = list(set(get_experiment_files(experiment, results_dir)))
    # all_scaleup_numbers.sort()
    # all_scaleup_numbers = [i for i in all_scaleup_numbers if i > 1]
    prefix = '{}/{}_'.format(results_dir, experiment)
    distr_speedup, compile_distr_speedup = collect_experiment_speedups(prefix, all_scaleup_numbers)

    output_diff = check_output_diff_correctness(prefix, all_scaleup_numbers)

    maximum_y = 0
    lines = []
    best_result = distr_speedup

    ## In the bigrams experiment we want to also have the alt_bigrams plot
    if(experiment == "bigrams" and not ("mini-split" in custom_scaleup_plots[experiment])):
        opt_prefix = '{}/{}_'.format(results_dir, "alt_bigrams")
        opt_distr_speedup, _ = collect_experiment_speedups(opt_prefix, all_scaleup_numbers)
        line, = ax.plot(all_scaleup_numbers, opt_distr_speedup, '-P', color='magenta', linewidth=0.5, label='Opt. Parallel')
        maximum_y = max(maximum_y, max(opt_distr_speedup))
        lines.append(line)

    if(not (experiment in custom_scaleup_plots) or
       (experiment in custom_scaleup_plots and
        "split" in custom_scaleup_plots[experiment])):
        try:
            auto_split_distr_speedup = collect_distr_experiment_speedup(prefix,
                                                                        'distr_auto_split.time',
                                                                        all_scaleup_numbers)
            line, = ax.plot(all_scaleup_numbers, auto_split_distr_speedup, '-D', color='tab:red', linewidth=0.5, label='Par + Split')
            lines.append(line)
            split_distr_speedup = collect_distr_experiment_speedup(prefix,
                                                                   'distr_split.time',
                                                                   all_scaleup_numbers)
            line, = ax.plot(all_scaleup_numbers, split_distr_speedup, '-h', color='brown', linewidth=0.5, label='Par + B. Split')
            lines.append(line)
            maximum_y = max(maximum_y, max(auto_split_distr_speedup))
            maximum_y = max(maximum_y, max(split_distr_speedup))
            best_result = split_distr_speedup
        except ValueError:
            pass

    if((experiment in custom_scaleup_plots and
        "mini-split" in custom_scaleup_plots[experiment])):
        try:
            auto_split_distr_speedup = collect_distr_experiment_speedup(prefix,
                                                                        'distr_auto_split.time',
                                                                        all_scaleup_numbers)
            line, = ax.plot(all_scaleup_numbers, auto_split_distr_speedup, '-^', linewidth=0.5)
            lines.append(line)
            maximum_y = max(maximum_y, max(auto_split_distr_speedup))
            best_result = auto_split_distr_speedup
        except ValueError:
            pass


    if(not (experiment in custom_scaleup_plots) or
       (experiment in custom_scaleup_plots and
        "eager" in custom_scaleup_plots[experiment])):
        line, = ax.plot(all_scaleup_numbers, distr_speedup, '-o', linewidth=0.5, label='Parallel')
        lines.append(line)
        maximum_y = max(maximum_y, max(distr_speedup))

    ## Add the no eager times if they exist
    if(not (experiment in custom_scaleup_plots) or
       (experiment in custom_scaleup_plots and
        "blocking-eager" in custom_scaleup_plots[experiment])):
        try:
            no_task_par_eager_distr_speedup = collect_distr_experiment_speedup(prefix,
                                                                               'distr_no_task_par_eager.time',
                                                                               all_scaleup_numbers)
            line, = ax.plot(all_scaleup_numbers, no_task_par_eager_distr_speedup, '-p', linewidth=0.5, label='Blocking Eager')
            lines.append(line)
            maximum_y = max(maximum_y, max(no_task_par_eager_distr_speedup))
        except ValueError:
            pass
    if(not (experiment in custom_scaleup_plots) or
       (experiment in custom_scaleup_plots and
        "no-eager" in custom_scaleup_plots[experiment])):
        try:
            no_eager_distr_speedup = collect_distr_experiment_speedup(prefix,
                                                                      'distr_no_eager.time',
                                                                      all_scaleup_numbers)
            line, = ax.plot(all_scaleup_numbers, no_eager_distr_speedup, '-^', linewidth=0.5, label='No Eager')
            lines.append(line)
            maximum_y = max(maximum_y, max(no_eager_distr_speedup))
            if((experiment in custom_scaleup_plots and
                not "eager" in custom_scaleup_plots[experiment])):
                best_result = no_eager_distr_speedup
        except ValueError:
            pass

    ## Set the ylim
    ax.set_ylim(top=maximum_y*1.15)

    return output_diff, lines, best_result


def collect_scaleup_times(experiment, results_dir):


    fig, ax = plt.subplots()

    ## Plot speedup
    ax.set_ylabel('Speedup')
    ax.set_xlabel('Level of Parallelism')

    # ax.plot(all_scaleup_numbers, total_distr_speedup, '-^', linewidth=0.5, label='+ Merge')
    # ax.plot(all_scaleup_numbers, all_scaleup_numbers, '-', color='tab:gray', linewidth=0.5, label='Ideal')
    all_scaleup_numbers = [2, 4, 8, 16, 32, 64]
    output_diff, _, _ = collect_scaleup_times_common(experiment, all_scaleup_numbers, results_dir, custom_scaleup_plots, ax)

    # plt.yscale("log")
    plt.xticks(all_scaleup_numbers[1:])
    plt.legend(loc='lower right')
    plt.title(pretty_names[experiment])


    plt.tight_layout()
    plt.savefig(os.path.join('../evaluation/plots', "{}_throughput_scaleup.pdf".format(experiment)),bbox_inches='tight')

    return output_diff

def plot_sort_with_baseline(results_dir):

    all_scaleup_numbers = [2, 4, 8, 16, 32, 64]
    sort_prefix = '{}/sort_'.format(results_dir)
    baseline_sort_prefix = '{}/baseline_sort/baseline_sort_'.format(results_dir)
    baseline_sort_opt_prefix = '{}/baseline_sort/baseline_sort_opt_'.format(results_dir)

    ## Collect all sort numbers
    seq_numbers, distr_numbers, _ = collect_experiment_scaleup_times(sort_prefix, all_scaleup_numbers)
    sort_distr_speedup = [safe_zero_div(seq_numbers[i], t) for i, t in enumerate(distr_numbers)]
    # sort_distr_speedup, _ = collect_experiment_speedups(sort_prefix, all_scaleup_numbers)
    baseline_sort_distr_speedup = collect_baseline_experiment_speedups(baseline_sort_prefix,
                                                                       [1] + [num*2
                                                                              for num in all_scaleup_numbers],
                                                                       seq_numbers[0])
    # baseline_sort_opt_distr_speedup = collect_baseline_experiment_speedups(baseline_sort_opt_prefix,
    #                                                                        [1] + all_scaleup_numbers[1:],
    #                                                                        seq_numbers[0])

    # output_diff = check_output_diff_correctness(prefix, all_scaleup_numbers)

    fig, ax = plt.subplots()

    ## Plot speedup
    ax.set_ylabel('Speedup')
    ax.set_xlabel('Level of Parallelism')
    ax.plot(all_scaleup_numbers, sort_distr_speedup, '-o', linewidth=0.5, label='Pash')
    ## Add the no eager times if they exist
    # try:
    #     no_task_par_eager_distr_speedup = collect_distr_experiment_speedup(sort_prefix,
    #                                                                         'distr_no_task_par_eager.time',
    #                                                                         all_scaleup_numbers)
    #     ax.plot(all_scaleup_numbers, no_task_par_eager_distr_speedup, '-p', linewidth=0.5, label='Pash - Blocking Eager')
    # except ValueError:
    #     pass

    try:
        no_eager_distr_speedup = collect_distr_experiment_speedup(sort_prefix,
                                                                   'distr_no_eager.time',
                                                                   all_scaleup_numbers)
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
    header += ['\\begin{tabular*}{\\textwidth}{l @{\\extracolsep{\\fill}} lllll}']
    header += ['\\toprule']
    if(not full):
        header += ['Script ~&~ Structure & Input &'
                   'Seq. Time & Highlights \\\\']
    else:
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

def generate_experiment_line(experiment, full=True):
    line = []
    line += [pretty_names[experiment], '~&~']
    line += [structures[experiment], '&']

    ## Collect and output the input size
    input_size = collect_format_input_size(experiment)
    line += [input_size, '&']

    ## Collect and output the sequential time for the experiment
    scaleup_numbers = [2, 16, 64]
    experiment_results_prefix = '{}/{}_'.format(RESULTS, experiment)
    seq_times, _, compile_times = collect_experiment_scaleup_times(experiment_results_prefix, scaleup_numbers)
    assert(len(seq_times) == 3)
    seq_time_seconds = format_time_seconds(seq_times[0])
    # seq_time_seconds = seq_times[0] / 1000
    line += [seq_time_seconds, '&']

    suffix='distr.time'
    if(experiment in ["spell", "bigrams", "double_sort"]):
        suffix='distr_auto_split.time'
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


def generate_tex_table(experiments):
    header = generate_table_header()
    lines = []
    for experiment in experiments:
        line = generate_experiment_line(experiment)
        # print(line)
        lines.append(line)
    data = "\n".join(lines)
    footer = generate_table_footer()
    table_tex = "\n".join([header, data, footer])
    tex_filename = os.path.join('../evaluation/plots', 'microbenchmarks-table.tex')
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


def collect_unix50_pipeline_scaleup_times(pipeline_number, unix50_results_dir, scaleup_numbers):
    prefix = '{}/unix50_pipeline_{}_'.format(unix50_results_dir, pipeline_number)
    distr_speedups, _ = collect_experiment_speedups(prefix, scaleup_numbers)
    absolute_seq_times, _, _ = collect_experiment_scaleup_times(prefix, scaleup_numbers)
    return (distr_speedups, absolute_seq_times)

def collect_unix50_pipeline_coarse_scaleup_times(pipeline_number, unix50_results_dir, scaleup_numbers):
    prefix = '{}/unix50_pipeline_{}_'.format(unix50_results_dir, pipeline_number)
    # distr_speedups, _ = collect_experiment_speedups(prefix, scaleup_numbers)
    no_eager_distr_speedup = collect_distr_experiment_speedup(prefix,
                                                              'distr_no_eager.time',
                                                              scaleup_numbers)
    absolute_seq_times = [read_total_time('{}{}_seq.time'.format(prefix, scaleup_numbers[0]))
                          for _ in scaleup_numbers]
    return (no_eager_distr_speedup, absolute_seq_times)


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

def make_unix50_bar_chart(all_results, scaleup_numbers, parallelism):

    ## TODO: Sort by speedup
    all_results.sort(key=lambda x: x[0][scaleup_numbers.index(parallelism)], reverse=True)
    ## Plot individual speedups
    individual_results = [distr_exec_speedup[scaleup_numbers.index(parallelism)]
                          for distr_exec_speedup, _ in all_results]
    absolute_seq_times_s = [absolute_seq[scaleup_numbers.index(parallelism)] / 1000
                            for _, absolute_seq in all_results]
    print("Unix50 individual speedups for {} parallelism:".format(parallelism), individual_results)
    mean = sum(individual_results) / len(individual_results)
    median = statistics.median(individual_results)
    geo_mean = math.exp(np.log(individual_results).sum() / len(individual_results))
    weighted_res = [i*a for i, a in zip(individual_results, absolute_seq_times_s)]
    weighted_avg = sum(weighted_res) / sum(absolute_seq_times_s)
    print("  Mean:", mean)
    print("  Median:", median)
    print("  Geometric Mean:", geo_mean)
    print("  Weighted Average:", weighted_avg)

    w = 0.2
    ind = np.arange(len(individual_results))
    speedup_color = 'tab:blue'
    seq_time_color = 'tab:red'

    fig, ax1 = plt.subplots()
    # print(fig.get_size_inches())
    fig.set_size_inches(12, 6)

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
    plt.savefig(os.path.join('../evaluation/plots', "unix50_individual_speedups_{}.pdf".format(parallelism)),bbox_inches='tight')

def make_coarse_unix50_bar_chart(all_results, scaleup_numbers, parallelism):
    ## Plot individual speedups
    individual_results = [distr_exec_speedup[scaleup_numbers.index(parallelism)]
                          for distr_exec_speedup, _ in all_results]
    absolute_seq_times_s = [absolute_seq[scaleup_numbers.index(parallelism)] / 1000
                            for _, absolute_seq in all_results]
    print("Unix50 individual speedups for {} parallelism:".format(parallelism), individual_results)
    mean = sum(individual_results) / len(individual_results)
    median = statistics.median(individual_results)
    geo_mean = math.exp(np.log(individual_results).sum() / len(individual_results))
    weighted_res = [i*a for i, a in zip(individual_results, absolute_seq_times_s)]
    weighted_avg = sum(weighted_res) / sum(absolute_seq_times_s)
    print("  Mean:", mean)
    print("  Median:", median)
    print("  Geometric Mean:", geo_mean)
    print("  Weighted Average:", weighted_avg)

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


def collect_unix50_scaleup_times(unix50_results_dir):
    files = [f for f in os.listdir(unix50_results_dir)]
    # print(files)
    pipeline_numbers = sorted(list(set([f.split('_')[2] for f in files])))
    # print(pipeline_numbers)

    scaleup_numbers = [2, 4, 8, 16]

    all_results = [collect_unix50_pipeline_scaleup_times(pipeline_number,
                                                         unix50_results_dir,
                                                         scaleup_numbers)
                   for pipeline_number in pipeline_numbers]
    # print(all_results)

    for parallelism in scaleup_numbers:
        make_unix50_bar_chart(all_results, scaleup_numbers, parallelism)


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
    plt.savefig(os.path.join('../evaluation/plots', "unix50_throughput_scaleup.pdf"),bbox_inches='tight')

def collect_unix50_coarse_scaleup_times(unix50_results_dir):
    files = [f for f in os.listdir(unix50_results_dir)]
    # print(files)
    pipeline_numbers = sorted(list(set([f.split('_')[2] for f in files])))
    # print(pipeline_numbers)

    scaleup_numbers = [16]

    all_results = [collect_unix50_pipeline_scaleup_times(pipeline_number,
                                                         unix50_results_dir,
                                                         scaleup_numbers)
                   for pipeline_number in pipeline_numbers]
    # print(all_results)

    for parallelism in scaleup_numbers:
        make_coarse_unix50_bar_chart(all_results, scaleup_numbers, parallelism)


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
    plt.savefig(os.path.join('../evaluation/plots', "unix50_coarse_throughput_scaleup.pdf"),bbox_inches='tight')


def plot_one_liners_tiling(results_dir):

    all_scaleup_numbers = [2, 4, 8, 16, 32, 64]
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

    confs = ["Par. + Split",
             "Par. + B. Split",
             "Parallel",
             "Blocking Eager",
             "No Eager",
             "Opt. Parallel"]

    fig = plt.figure()
    gs = fig.add_gridspec(2, 5, hspace=0.05)
    # fig.suptitle('')

    total_lines = []
    averages = [[] for _ in all_scaleup_numbers]
    ## Plot microbenchmarks
    for i, experiment in enumerate(experiments):
        ax = fig.add_subplot(gs[i])
        _, lines, best_result = collect_scaleup_times_common(experiment, all_scaleup_numbers, results_dir, custom_scaleup_plots, ax)
        if(experiment == "double_sort"):
            total_lines = lines + total_lines
        elif(experiment == "bigrams"):
            total_lines += [lines[0]]
        ax.set_xticks(all_scaleup_numbers[1:])
        ax.text(.5,.91,pretty_names[experiment],
        horizontalalignment='center',
        transform=ax.transAxes)
        # ax.set_yticks([])
        fig.add_subplot(ax)

        ## Update averages
        for i, res in enumerate(best_result):
            averages[i].append(res)


    axs = fig.get_axes()
    for ax in axs:
        if(ax.is_first_col()):
            ax.set_ylabel('Speedup')
        if(ax.is_last_row()):
            ax.set_xlabel('Level of Parallelism')
        if(ax.is_first_row()):
            ax.set_xticklabels([])
        # ax.label_outer()

    plt.legend(total_lines, confs, loc='lower right', fontsize=16)
    # plt.title(pretty_names[experiment])

    fig.set_size_inches(27, 8.2)
    plt.tight_layout()
    plt.savefig(os.path.join('../evaluation/plots', "tiling_throughput_scaleup.pdf"),bbox_inches='tight')

    ## Print average, geo-mean
    one_liner_averages = [sum(res)/len(res) for res in averages]
    geo_means = [math.exp(np.log(res).sum() / len(res))
                 for res in averages]
    print("One-liners Aggregated results:")
    print(" |-- Averages:", one_liner_averages)
    print(" |-- Geometric Means:", geo_means)

def plot_less_one_liners_tiling(results_dir):

    all_scaleup_numbers = [2, 4, 8, 16, 32, 64]
    experiments = ["minimal_grep",
                   "minimal_sort",
                   "topn",
                   "wf",
                   "spell",
                   "diff",
                   "bigrams",
                   "set-diff",
                   "shortest_scripts"]

    # confs = ["Parallel"]

    fig = plt.figure()
    gs = fig.add_gridspec(3, 3, hspace=0.05)
    # fig.suptitle('')

    total_lines = []
    averages = [[] for _ in all_scaleup_numbers]
    ## Plot microbenchmarks
    for i, experiment in enumerate(experiments):
        ax = fig.add_subplot(gs[i])
        _, lines, best_result = collect_scaleup_times_common(experiment, all_scaleup_numbers, results_dir, coarse_custom_scaleup_plots, ax)
        if(experiment == "double_sort"):
            total_lines = lines + total_lines
        elif(experiment == "bigrams"):
            total_lines += [lines[0]]
        ax.set_xticks(all_scaleup_numbers[1:])
        ax.text(.5,.91,pretty_names[experiment],
        horizontalalignment='center',
        transform=ax.transAxes)
        # ax.set_yticks([])
        fig.add_subplot(ax)

        ## Update averages
        for i, res in enumerate(best_result):
            averages[i].append(res)


    axs = fig.get_axes()
    for ax in axs:
        if(ax.is_first_col()):
            ax.set_ylabel('Speedup')
        if(ax.is_last_row()):
            ax.set_xlabel('Level of Parallelism')
        if(not ax.is_last_row()):
            ax.set_xticklabels([])
        # ax.label_outer()

    # plt.legend(total_lines, confs, loc='lower right', fontsize=16)
    # plt.title(pretty_names[experiment])

    fig.set_size_inches(22.5, 12.5)
    plt.tight_layout()
    plt.savefig(os.path.join('../evaluation/plots', "coarse_tiling_throughput_scaleup.pdf"),bbox_inches='tight')

    ## Print average, geo-mean
    one_liner_averages = [sum(res)/len(res) for res in averages]
    geo_means = [math.exp(np.log(res).sum() / len(res))
                 for res in averages]
    print("Coarse One-liners Aggregated results:")
    print(" |-- Averages:", one_liner_averages)
    print(" |-- Geometric Means:", geo_means)



def format_wrong_output(output_diff, experiment, mode):
    if(len(output_diff) > 0):
        formatted_output_diff = 'Output for {} -- {} {} is wrong'.format(experiment, output_diff, mode)
        return [' !! -- WARNING -- !! {}'.format(formatted_output_diff)]
    else:
        return []

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
for experiment in experiments:
    output_diff = collect_scaleup_times(experiment, RESULTS)
    diff_results += format_wrong_output(output_diff[0], experiment, "parallel")
    diff_results += format_wrong_output(output_diff[1], experiment, "parallel no-eager")
    diff_results += format_wrong_output(output_diff[2], experiment, "parallel no-task-par-eager")

plot_one_liners_tiling(RESULTS)

plot_less_one_liners_tiling(RESULTS)

## Generate Tex table for microbenchmarks
generate_tex_table(experiments)

generate_tex_coarse_table(coarse_experiments)


## Plot Unix50
collect_unix50_scaleup_times(UNIX50_RESULTS)
collect_unix50_coarse_scaleup_times(UNIX50_RESULTS)


## Plot sort against sort with parallel-flag
plot_sort_with_baseline(RESULTS)

print("\n".join(diff_results))
