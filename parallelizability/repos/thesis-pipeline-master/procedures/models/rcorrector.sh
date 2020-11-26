#! /usr/bin/bash

# Exit on First Error - to Prevent Invalid File Modifications
# set -o errexit

# Load ~/.bash_profile if Not Found
if [ -z $PIPELINE_HOME ]; then
    echo "Reloading ~/.bash_profile"
    source ~/.bash_profile
fi

#
# Assign Arguments
# 

for i in "$@"
    do case $i in

    # Standard Arguments

        # Access & Write Files With This Prefix
        -f=*|--fileprefix=*)
        fileprefix="${i#*=}"
        shift
        ;;

        # Access & Write Files With This Subset
        -s=*|--subset=*)
        subset="${i#*=}"
        shift
        ;;

        # Access & Write Files With This Condition
        -c=*|--condition=*)
        condition="${i#*=}"
        shift
        ;;

        # Access & Write Files With This Experiment
        -x=*|--experiment=*)
        experiment="${i#*=}"
        shift
        ;;

        # Access & Write Files With This Parameter Set
        -p=*|--parameters=*)
        parameters="${i#*=}"
        shift
        ;;

        # Access & Write Files With This Read Group
        -g=*|--readgroup=*)
        readgroup="${i#*=}"
        shift 
        ;;

    # Optional Arguments With Defaults

        # Number of Cores to Use
        -n=*|--ncores=*)
        ncoresOpt="${i#*=}"
        shift
        ;;

        # Per Core Memory Requirement
        -m=*|--memory=*)
        memoryOpt="${i#*=}"
        shift
        ;;

    # Invalid Argument Handler

        *)
        # invalid option
        printf "Invalid/Unused Parameter: $i"
        ;;
        
    esac
done

# Get Max Allowable Memory
allocMemory=${memory//[GgMmKk]/}
allocSize=${memory//[0-9]/}
allocMax=$((allocMemory * ncores))
maxMemory=$((allocMemory * ncores))$allocSize

# Set Directories
dataDir=$PIPELINE_HOME/$subset
paramDir=$PIPELINE_HOME/$subset/model/$experiment/param/$parameters

# 
# Run Rcorrector
# 

# State Check - Run Block if it Has Not Already Been Executed Successfully
state="$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup:RCORRECTOR:1"
if !(has_state $state); then

    if [ "$parameters" = "default" ]; then
        
        #
        # Default Parameters
        #

        format_status "Running Rcorrector - Default Parameters"
        # Define Command
        call="rcorrector \
        -i $dataDir/fastq/split/$fileprefix.$subset.$condition.$readgroup.fastq \
        -k 31 \
        -t $ncores \
        -od $paramDir/modeled"
        # Print Command & Call
        format_status "Command:\n$call"
        eval $call

    elif [ "$parameters" = "custom" ]; then

        #
        # Custom Parameters
        #

        format_status "Running Rcorrector - Custom Parameters"
        # Define Command
        call="rcorrector \
        -i $dataDir/fastq/split/$fileprefix.$subset.$condition.$readgroup.fastq \
        -k 31 \
        -t $ncores \
        -od $paramDir/modeled"
        # Print & Call
        format_status "Command:\n$call"
        eval $call

    fi
    
    # Update State on Exit
    status=$?
    put_state $status $state
    format_status "Rcorrector ($parameters $readgroup) Complete"
    return $status

fi

# Usage: perl ./run_rcorrector.pl [OPTIONS]
# OPTIONS:
# Required parameters:
# 	-s seq_files: comma separated files for single-end data sets
# 	-1 seq_files_left: comma separated files for the first mate in the paried-end data sets
# 	-2 seq_files_right: comma separated files for the second mate in the paired-end data sets
# 	-i seq_files_interleaved: comma sperated files for interleaved paired-end data sets
# Other parameters:
# 	-k kmer_length (<=32, default: 23)
# 	-od output_file_directory (default: ./)
# 	-t number_of_threads (default: 1)
# 	-maxcorK INT: the maximum number of correction within k-bp window (default: 4)
# 	-wk FLOAT: the proportion of kmers that are used to estimate weak kmer count threshold, lower for more divergent genome (default: 0.95)
# 	-ek expected_number_of_kmers: does not affect the correctness of program but affect the memory usage (default: 100000000)
# 	-stdout: output the corrected reads to stdout (default: not used)
# 	-verbose: output some correction information to stdout (default: not used)
# 	-stage INT: start from which stage (default: 0)
# 		0-start from begining(storing kmers in bloom filter);
# 		1-start from count kmers showed up in bloom filter;
# 		2-start from dumping kmer counts into a jf_dump file;
# 		3-start from error correction.