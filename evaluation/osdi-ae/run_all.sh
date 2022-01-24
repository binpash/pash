#!/bin/bash
# go to benchmark directory
cd ../benchmarks
setup_flags='--small'
echo 'Downloading Package Deps'
# install the benchmark deps
bash web-index/input/install-deps.sh 1> /dev/null 2> /dev/null
bash dependency_untangling/input/install-deps.sh 1> /dev/null 2> /dev/null
echo 'Deps downloaded'

run_bash() {
    ## This script is necessary to ensure that sourcing happens with bash
    source run.seq.sh
    bench_len=$((${#PASH_BENCHMARK[@]} -1))
    array_len=$((${#PASH_ALL_FLAGS[@]} -1))
    for i in $(seq 0 $bench_len)
    do
        export IN=""
        export IN_PRE=""
        bench=${PASH_BENCHMARK[$i]}
        cd $bench
        # remove all the time files
        rm -rf outputs output pash_logs par.res
        cd ..
        echo 'Running bash:' ${bench}
        bdir=${OUT}/bash/${bench}
        mkdir -p ${bdir}
        # run the benchmark
        ${bench} ${setup_flags}
        # copy the output
        cp -r ${bench}/outputs/* ${bdir}/
        # copy the time file
        cp ${bench}/seq.res ${bdir}/
        cd ${bench}
        rm -rf outputs output pash_logs *.res
        cd ..
    done
}

run_bench() {
    ## This script is necessary to ensure that sourcing happens with bash
    source run.par.sh
    bench_len=$((${#PASH_BENCHMARK[@]} -1))
    array_len=$((${#PASH_ALL_FLAGS[@]} -1))
    for i in $(seq 0 $bench_len)
    do
        export IN=""
        export IN_PRE=""
        bench=${PASH_BENCHMARK[$i]}
        cd $bench
        # remove all the time files
        rm -rf outputs output pash_logs par.res
        cd ..
        for j in $(seq 0 $array_len)
        do
            export mode=${PASH_MODE[$j]}
            export PASH_FLAGS=${PASH_ALL_FLAGS[$j]}
            ${bench}_pash ${setup_flags}
            pdir=${OUT}/${mode}/${bench}
            mkdir -p ${pdir}
            # move the folder to our dest
            mv ${bench}/outputs ${pdir}/
            # copy the time file
            mv ${bench}/par.res ${pdir}/
        done
        cd ${bench}
        rm -rf outputs output pash_logs *.res
        cd ..
    done
}

OUT=artifact/
rm -rf ${OUT}
mkdir -p ${OUT}

export PASH_ALL_FLAGS=(" "
                       "--r_split --dgsh_tee --r_split_batch_size 1000000 --parallel_pipelines --profile_driven")
export PASH_BENCHMARK=("oneliners" "unix50" "analytics-mts" "poets" "max-temp" "web-index" "dependency_untangling")
export PASH_MODE=("pash" 
                  "blish")

echo 'Running all benchmark for bash'
time run_bash

echo 'Running Blish benchmarks'
time run_bench 

##### Figure 6
export PASH_ALL_FLAGS=("--r_split --dgsh_tee --r_split_batch_size 1000000" 
                       "--r_split --dgsh_tee --r_split_batch_size 1000000 --parallel_pipelines" )
export PASH_BENCHMARK=("poets" "max-temp" "dependency_untangling")
export PASH_MODE=("blish_no_prof_no_du" 
                  "blish_no_prof")

time run_bench 

##### Figure 7
export PASH_ALL_FLAGS=(
"--dgsh_tee --parallel_pipelines --profile_driven" )
export PASH_BENCHMARK=("oneliners" "unix50" "analytics-mts" "max-temp" "web-index")
export PASH_MODE=("blish_no_comm")

time run_bench 

# kill the hanging processes 
pkill -f cat
