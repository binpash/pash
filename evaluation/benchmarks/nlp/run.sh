#!/bin/bash

export SUITE_DIR=$(realpath $(dirname "$0"))
export TIMEFORMAT=%R
cd $SUITE_DIR

if [[ "$1" == "--small" ]]; then
    echo "Using small input"
    export ENTRIES=10
    export IN="$SUITE_DIR/inputs/pg-small"
else
    echo "Using default input"
    export ENTRIES=120
    export IN="$SUITE_DIR/inputs/pg"
fi

names_scripts=(
    "1syllable_words;6_4"
    "2syllable_words;6_5"
    "4letter_words;6_2"
    "bigrams_appear_twice;8.2_2"
    "bigrams;4_3"
    "compare_exodus_genesis;8.3_3"
    "count_consonant_seq;7_2"
    "count_morphs;7_1"
    "count_trigrams;4_3b"
    "count_vowel_seq;2_2"
    "count_words;1_1"
    "find_anagrams;8.3_2"
    "merge_upper;2_1"
    "sort;3_1"
    "sort_words_by_folding;3_2"
    "sort_words_by_num_of_syllables;8_1"
    "sort_words_by_rhyming;3_3"
    "trigram_rec;6_1" # was initially commented out
    "uppercase_by_token;6_1_1"
    "uppercase_by_type;6_1_2"
    "verses_2om_3om_2instances;6_7"
    "vowel_sequencies_gr_1K;8.2_1"
    "words_no_vowels;6_3"
)

mkdir -p "outputs"
all_res_file="./outputs/nlp.res"
> $all_res_file

# time_file stores the time taken for each script
# mode_res_file stores the time taken and the script name for every script in a mode (e.g. bash, pash, dish, fish)
# all_res_file stores the time taken for each script for every script run, making it easy to copy and paste into the spreadsheet
nlp() {
    mkdir -p "outputs/$1"
    mode_res_file="./outputs/$1/nlp.res"
    > $mode_res_file

    echo executing nlp $1 $(date) | tee -a $mode_res_file $all_res_file

    for name_script in ${names_scripts[@]}
    do
        IFS=";" read -r -a name_script_parsed <<< "${name_script}"
        name="${name_script_parsed[0]}"
        script="${name_script_parsed[1]}"
        script_file="./scripts/$script.sh"
        output_dir="./outputs/$1/$script/"
        output_file="./outputs/$1/$script.out"
        time_file="./outputs/$1/$script.time"
        log_file="./outputs/$1/$script.log"
        # output_file contains "done" when run successfully. The real outputs are under output_dir/
        if [[ "$1" == "bash" ]]; then
            (time bash $script_file $output_dir > $output_file) 2> $time_file
        fi

        cat "${time_file}" >> $all_res_file
        echo "$script_file $(cat "$time_file")" | tee -a $mode_res_file
    done
}

nlp "bash"
