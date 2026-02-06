#!/bin/bash

# Input matrix for evaluation/benchmarks/unix50.
# Mirrors the original unix50 input selections (including commented candidates).
set_leash_benchmark_inputs() {
    if [[ "$*" == *"--small"* ]]; then
        SCRIPT_INPUT_WIDTH=(
            # "1.sh:1_1M.txt"
            # # "2.sh:1_1M.txt"
            # "3.sh:1_1M.txt"
            # "4.sh:1_1M.txt"
            # "5.sh:2_1M.txt"
            "6.sh:3_1M.txt"
            # "7.sh:4_1M.txt"
            # "8.sh:4_1M.txt"
            # "9.sh:4_1M.txt"
            # "10.sh:4_1M.txt"
            # # "11.sh:4_1M.txt"
            # "13.sh:5_1M.txt"
            # # "14.sh:6_1M.txt"
            # "15.sh:7_1M.txt"
            # "17.sh:7_1M.txt"
            # "18.sh:8_1M.txt"
            # "19.sh:8_1M.txt"
            # "20.sh:8_1M.txt"
            # # "21.sh:8_1M.txt"
            # "23.sh:9.1_1M.txt"
            # "24.sh:9.2_1M.txt"
            # "25.sh:9.3_1M.txt"
            # "26.sh:9.4_1M.txt"
            # "28.sh:9.6_1M.txt"
            # "29.sh:9.7_1M.txt"
            # "30.sh:9.8_1M.txt"
            # "31.sh:9.9_1M.txt"
            # "32.sh:10_1M.txt"
            # "33.sh:10_1M.txt"
            # "35.sh:11_1M.txt"
        )
    elif [[ "$*" == *"--medium"* ]]; then
        SCRIPT_INPUT_WIDTH=(
            # "1.sh:1_1M.txt"
            # # "2.sh:1_1M.txt"
            # "3.sh:1_1M.txt"
            # "4.sh:1_1M.txt"
            # "5.sh:2_1M.txt"
            "6.sh:3_5G.txt"
            # "7.sh:4_1M.txt"
            # "8.sh:4_1M.txt"
            # "9.sh:4_1M.txt"
            # "10.sh:4_1M.txt"
            # # "11.sh:4_1M.txt"
            # "13.sh:5_1M.txt"
            # # "14.sh:6_1M.txt"
            # "15.sh:7_1M.txt"
            # "17.sh:7_1M.txt"
            # "18.sh:8_1M.txt"
            # "19.sh:8_1M.txt"
            # "20.sh:8_1M.txt"
            # # "21.sh:8_1M.txt"
            # "23.sh:9.1_1M.txt"
            # "24.sh:9.2_1M.txt"
            # "25.sh:9.3_1M.txt"
            # "26.sh:9.4_1M.txt"
            # "28.sh:9.6_1M.txt"
            # "29.sh:9.7_1M.txt"
            # "30.sh:9.8_1M.txt"
            # "31.sh:9.9_1M.txt"
            # "32.sh:10_1M.txt"
            # "33.sh:10_1M.txt"
            # "35.sh:11_1M.txt"
        )
    elif [[ "$*" == *"--large"* ]]; then
        SCRIPT_INPUT_WIDTH=(
            # "1.sh:1_20G.txt"
            # "3.sh:1_20G.txt" # Error
            # "2.sh:1_20G.txt" # 4096 Mem usage
            # "4.sh:1_20G.txt" # 4096 Mem usage
            # "5.sh:2_20G.txt"
            # "6.sh:3_20G.txt"
            # "7.sh:4_20G.txt"
            # "8.sh:4_20G.txt"
            # "9.sh:4_20G.txt"
            "10.sh:4_20G.txt"
            # "11.sh:4_20G.txt" # sort
            # 12.sh:4_20G.txt:16
            # "13.sh:5_20G.txt"
            # "14.sh:6_20G.txt" # sort
            # "15.sh:7_20G.txt"
            # "17.sh:7_20G.txt" # sort
            # "18.sh:8_20G.txt"
            "19.sh:8_20G.txt"
            # # "20.sh:8_20G.txt"
            # "21.sh:8_20G.txt"
            # "23.sh:9.1_20G.txt:32"
            # "24.sh:9.2_20G.txt"
            # "25.sh:9.3_20G.txt"
            # "26.sh:9.4_20G.txt"
            # "28.sh:9.6_20G.txt:16"
            # "29.sh:9.7_20G.txt"
            # "30.sh:9.8_20G.txt"
            # "31.sh:9.9_20G.txt"
            # "32.sh:10_20G.txt"
            # "33.sh:10_20G.txt"
            # "35.sh:11_20G.txt"
            # "34.sh:10_20G.txt"
        )
    else
        SCRIPT_INPUT_WIDTH=(
            # "1.sh:1_1G.txt"
            # "3.sh:1_1G.txt" # Error
            # "2.sh:1_1G.txt" # 4096 Mem usage
            # "4.sh:1_1G.txt" # 4096 Mem usage
            # "5.sh:2_1G.txt"
            "6.sh:3_1G.txt"
            # "7.sh:4_1G.txt"
            # "8.sh:4_1G.txt"
            # "9.sh:4_1G.txt"
            # "10.sh:4_1G.txt" # sort
            # "11.sh:4_1G.txt" # sort
            # 12.sh:4_1G.txt:16
            # "13.sh:5_1G.txt"
            # "14.sh:6_1G.txt" # sort
            # "15.sh:7_1G.txt"
            # "17.sh:7_1G.txt" # sort
            # "18.sh:8_1G.txt"
            # "19.sh:8_1G.txt"
            # # "20.sh:8_1G.txt"
            # "21.sh:8_1G.txt"
            # "23.sh:9.1_1G.txt:32"
            # "24.sh:9.2_1G.txt"
            # "25.sh:9.3_1G.txt"
            # "26.sh:9.4_1G.txt"
            # "28.sh:9.6_1G.txt:16"
            # "29.sh:9.7_1G.txt"
            # "30.sh:9.8_1G.txt"
            # "31.sh:9.9_1G.txt"
            # "32.sh:10_1G.txt"
            # "33.sh:10_1G.txt"
            # "35.sh:11_1G.txt"
            # "34.sh:10_1G.txt"
        )
    fi
}
