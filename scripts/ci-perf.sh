#! /usr/bin/env bash

# Run performance tests

main() {
    set -Eex;

    local pash_d="$(get_pash_dir)";

    cd "$pash_d";
    git fetch;
    local initial_revision="$(get_revision HEAD)";
    local latest_main_revision="$(get_revision main)";
    local revision="${1:-$latest_main_revision}";

    local output_dir="${2:-/tmp/results}";
    local output_revision_directory="${output_dir}/$revision";
    echo "Will write to $output_revision_directory";

    # For reproducibility.
    trap "git checkout '$initial_revision'" EXIT

    # Use subshell for new working directory and
    # visual distinction in `set -e`
    echo "Running performance tests for $revision"
    (git checkout "$revision" && \
     build_pash_runtime && \
     run_performance_test_suites);

    mkdir -p "$output_revision_directory";
    cp -r "$pash_d/evaluation/results/." "$output_revision_directory/"

    # The code to build the summary file might not be in the commit
    # used to run the tests.
    git checkout "$latest_main_revision";

    echo "Summarizing results";
    local eurosys_tests='bigrams,diff,minimal_grep,minimal_sort,set-diff,spell,topn,wf'
    summarize_perf_suite "EuroSys One-liners" \
                         "$revision" \
                         "${output_revision_directory}/eurosys_small" \
                         "$eurosys_tests" \
                         "2" \
                         "distr_auto_split" \
                         "${output_dir}/summary_eurosys_small"

    # Generate index page so others can review available summaries
    # through web server.
    cd "${output_dir}"
    ls summary_* > index;
    cd -
}


build_pash_runtime() {
    make -C "$(get_pash_dir)/runtime";
}

get_pash_dir() {
    local here="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )";
    git -C "$here" rev-parse --show-toplevel;
}

get_revision() {
    git rev-parse --short "${1:-HEAD}";
}

run_performance_test_suites() {
    local pash_d=$(get_pash_dir);
    cd "$pash_d/evaluation/eurosys";
    ./execute_eurosys_one_liners.sh -s
    # ./execute_unix_benchmarks.sh -l
    # ./execute_baseline_sort.sh
    # ./execute_max_temp_dish_evaluation.sh
    # ./execute_web_index_dish_evaluation.sh
}

summarize_perf_suite() {
    local heading="$1";
    local revision="$2";
    local input_dir="$3";
    local tests="$4";
    local width="$5";
    local variant="$6";
    local summary_file="$7";
    local cell_fmt='%-20s';

    IFS=',' read -ra test_array <<< "$tests";

    # When starting a summary file, include a header.
    if [[ ! -f "$summary_file" ]]; then
        (
            printf "$heading (width=$width variant=$variant)\n";
            printf "$cell_fmt" 'revision';
            for t in "${test_array[@]}"; do
                printf "$cell_fmt" "$t";
            done;
            printf '\n';
        ) > "$summary_file";
    fi

    # Add a row of test data.
    printf "$cell_fmt" "$revision" >> "$summary_file";
    for t in "${test_array[@]}"; do
        local perf_file="${input_dir}/${t}_${width}_${variant}.time";
        echo "Summarizing $perf_file";
        printf "$cell_fmt" $(summarize_perf_file "$perf_file") >> "$summary_file";
    done
    printf '\n' >> "$summary_file";
}

print_pash_execution_time() {
    LC_NUMERIC='C' \
        cat "$1" | \
        grep 'Execution time: ' | \
        sed 's/[^0-9\.]//g' | \
        awk '{s+=sprintf("%f", $1)}END{printf "%.4f",s}';
}

print_user_time() {
    local time_string="$(egrep 'user[^m]+m[0-9\.]+s' "$1" | sed 's/^[^0-9]+//g')";
    local seconds="$(echo "$time_string" | sed -nr 's/.*m([^s]+)s/\1/p')"
    local minutes="$(echo "$time_string" | sed -nr 's/^[^0-9]+([0-9\.]+)m.*/\1/p')";
    echo "scale=4; ($minutes * 60) + $seconds" | bc;
}

summarize_perf_file() {
    local perf_file="$1";
    read -a data < <(split_perf_file_name "$perf_file");

    local test="${data[0]}";
    local width="${data[1]}";
    local variant="${data[2]}";

    if [[ "$variant" == 'seq' ]]; then
        printf "%ss" "$(print_user_time "$1")";
    elif [[ -f "$(make_perf_file_name "$test" "$width" "seq")" ]]; then
        local ptime="$(print_pash_execution_time "$1")";
        local utime="$(print_user_time "$1")";
        printf "%ss,x%s" "$ptime" "$(echo "scale=4; $utime / $ptime" | bc)";
    else
        print_pash_execution_time "$1";
    fi
}

split_perf_file_name() {
    if [[ "$(basename $1)" =~ (.*)_([0-9]+)_(.*).time$ ]]; then
        echo "${BASH_REMATCH[@]:1}";
        return 0
    else
        return 1
    fi
}

make_perf_file_name() {
    local name="$1";
    local width="$2";
    local variant="$3";
    echo "${name}_${width}_${variant}.time";
}


(return 0 2>/dev/null) || main "$@"
