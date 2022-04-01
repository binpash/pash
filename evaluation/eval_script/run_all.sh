#!/bin/bash
RES_FOLDER=${PWD}/eval_results/run
# go to benchmark directory
cd ${PASH_TOP}/evaluation/benchmarks
# use the small input for the benchmarks
setup_flags='--small'

# run all the scripts using bash
run_bash() {
    ## This script is necessary to ensure that sourcing happens with bash
    source run.seq.sh
    bench_len=$((${#PASH_BENCHMARK[@]} -1))
    array_len=$((${#PASH_ALL_FLAGS[@]} -1))
    for i in $(seq 0 $bench_len)
    do
        export IN=
        export IN_PRE=
        bench=${PASH_BENCHMARK[$i]}
        echo 'Running bash:' ${bench}
        bdir=${RES_FOLDER}/bash/${bench}
        mkdir -p ${bdir}
        # run the benchmark
        ${bench} ${setup_flags}
        # copy the time file
        mv ${bench}/seq.res ${bdir}/
    done
}

# run all the scripts using different configurations of PaSh JIT/PaSh AOT
run_bench() {
    ## This script is necessary to ensure that sourcing happens with bash
    source run.par.sh
    bench_len=$((${#PASH_BENCHMARK[@]} -1))
    array_len=$((${#PASH_ALL_FLAGS[@]} -1))
    for i in $(seq 0 $bench_len)
    do
        bench=${PASH_BENCHMARK[$i]}
        # remove all the time files
        for j in $(seq 0 $array_len)
        do
            export IN=
            export IN_PRE=
            export mode=${PASH_MODE[$j]}
            export PASH_FLAGS=${PASH_ALL_FLAGS[$j]}
            pdir=${RES_FOLDER}/${mode}/${bench}
            ${bench}_pash ${setup_flags}
            mkdir -p ${pdir}
            # move the folder to our dest
            rm -rf ${bench}/outputs
            # copy the time file
            mv ${bench}/par.res ${pdir}/
        done
    done
}

function run_all_benchmarks() {
  # generate output folder for each run
  export RES_FOLDER=$1
  # clean previous runs
  rm -rf ${RES_FOLDER}
  mkdir -p ${RES_FOLDER}
  cd ${PASH_TOP}/evaluation/benchmarks
  # remove all res files from previous runs
  find . -type d -name "outputs" | xargs rm -rf
  # do not remove any input from the node_modules dataset
  find . -type d -not -path "*/node_modules/*" -name "output" | xargs rm -rf 
  find . -type d -name "pash_logs" | xargs rm -rf
  find . -type f -name "*.res" | xargs rm -f
  # start preparing from execution
  export PASH_ALL_FLAGS=(" "
                         "--r_split --dgsh_tee --r_split_batch_size 1000000 --parallel_pipelines --profile_driven")
  export PASH_BENCHMARK=("oneliners" "unix50" "analytics-mts" "nlp" "max-temp" "web-index" "dependency_untangling")
  export PASH_MODE=("pash_aot" 
                    "pash_jit")

  echo 'Running all benchmark for bash'
  time run_bash
  echo 'Running PaSh JIT/PaSh AOT benchmarks'
  time run_bench 

  ##### Figure 6
  export PASH_ALL_FLAGS=("--r_split --dgsh_tee --r_split_batch_size 1000000" 
                         "--r_split --dgsh_tee --r_split_batch_size 1000000 --parallel_pipelines" )
  export PASH_BENCHMARK=("nlp" "max-temp" "dependency_untangling")
  export PASH_MODE=("pash_jit_no_prof_no_du" 
                    "pash_jit_no_prof")

  time run_bench 

  ##### Figure 7
  export PASH_ALL_FLAGS=(
  #"--dgsh_tee  # omitted until it's fixed
  "--parallel_pipelines --profile_driven" )
  export PASH_BENCHMARK=("oneliners" "unix50" "analytics-mts" "max-temp" "web-index")
  export PASH_MODE=("pash_jit_no_comm")

  time run_bench 

  # kill the hanging processes 
  pkill -f cat
}
out_folder=${RES_FOLDER}
# run all the tests three times
# the results will be stored at: $PASH_TOP/evaluation/eval_results/run{1,2,3}
run_all_benchmarks ${out_folder}1
run_all_benchmarks ${out_folder}2
run_all_benchmarks ${out_folder}3
