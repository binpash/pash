#!/bin/bash

# time: print real in seconds, to simplify parsing
TIMEFORMAT="%3R" # %3U %3S"

if [[ -z "$PASH_TOP" ]]; then
  echo "Need to provide PASH_TOP, possibly $(git rev-parse --show-toplevel)" 1>&2
  exit 1
fi

eval_dir="$PASH_TOP/evaluation/benchmarks/runtime-overhead/"

bash_outputs_suffix="bash.out"
par_outputs_suffix="par.out"
outputs_dir="$eval_dir/outputs"
pash_logs_dir="$eval_dir/pash_logs"

mkdir -p "$outputs_dir"
mkdir -p "$pash_logs_dir"

times_file="$eval_dir/time.res"

script_name="for-echo"
script="${script_name}.sh"

# The number of loop iterations
export N=100

printf -v pad %40s

## Bash
bash_outputs_file="${outputs_dir}/${script_name}.${bash_outputs_suffix}"
config="bash:${pad}"
config=${config:0:40}
echo "${config}" $({ time bash ${script} > "$bash_outputs_file" ; } 2>&1) | tee "$times_file"

run_pash()
{
    local config="$1"
    local PASH_FLAGS="$2"
    config_padded="$config:${pad}"
    config_padded=${config_padded:0:40}
    par_outputs_file="${outputs_dir}/${script_name}.${config}.${par_outputs_suffix}"
    pash_log="${pash_logs_dir}/${script_name}.${config}.pash.log"

    ## We don't want -d 1 since it adds overhead!
    echo "${config_padded}" $({ time "$PASH_TOP/pa.sh" $PASH_FLAGS  --log_file "${pash_log}" ${script} > "$par_outputs_file"; } 2>&1) | tee -a "$times_file"
    diff -q "$bash_outputs_file" "$par_outputs_file"
}

config="PaSh_no_daemon"
PASH_FLAGS="--no_daemon"

run_pash "$config" "$PASH_FLAGS"

config="PaSh_daemon_bash_mirror"
PASH_FLAGS="--expand_using_bash_mirror"

run_pash "$config" "$PASH_FLAGS"

config="PaSh_daemon"
PASH_FLAGS=""

run_pash "$config" "$PASH_FLAGS"

config="PaSh_daemon_fifos"
PASH_FLAGS="--daemon_communicates_through_unix_pipes"

run_pash "$config" "$PASH_FLAGS"

config="PaSh_daemon_par_pipelines"
PASH_FLAGS="--parallel_pipelines"

run_pash "$config" "$PASH_FLAGS"

config="PaSh_daemon_par_pipelines_fifos"
PASH_FLAGS="--parallel_pipelines --daemon_communicates_through_unix_pipes"

run_pash "$config" "$PASH_FLAGS"
