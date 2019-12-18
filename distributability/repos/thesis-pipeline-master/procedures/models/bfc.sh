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
# Run BFC
# 

# State Check - Run Block if it Has Not Already Been Executed Successfully
state="$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup:BFC:1"
if !(has_state $state); then

    # No Test for Paired / Single Reads Required

    if [ "$parameters" = "default" ]; then
        
        #
        # Default Parameters
        #

        format_status "Running BFC - Default Parameters"
        # Call Error Model - Size Parameter = Approx. 3 Gigabase Human Genome (hg19)
        # Define Command
        call="$BFC \
        -s 3g \
        -t $ncores \
        $dataDir/fastq/split/$fileprefix.$subset.$condition.$readgroup.fastq \
        > $paramDir/modeled/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.fastq"
        # Print & Call
        format_status "Command:\n$call"
        eval $call


    elif [ "$parameters" = "custom" ]; then

        #
        # Custom Parameters
        #

        format_status "Running BFC - Custom Parameters"
        # Call Error Model - Size Parameter = Approx. 3 Gigabase Human Genome (hg19)
        # Define Command
        call="$BFC \
        -s 3g \
        -t $ncores \
        $dataDir/fastq/split/$fileprefix.$subset.$condition.$readgroup.fastq \
        > $paramDir/modeled/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.fastq"
        # Print & Call
        format_status "Command:\n$call"
        eval $call

    fi
    
    # Update State on Exit
    status=$?
    put_state $status $state
    format_status "BFC ($parameters $readgroup) Complete"
    return $status

fi

# ./bfc
# Usage: bfc [options] <to-count.fq> [to-correct.fq]
# Options:
#   -s FLOAT     approx genome size (k/m/g allowed; change -k and -b) [unset]
#   -k INT       k-mer length [33]
#   -t INT       number of threads [1]
#   -b INT       set Bloom filter size to pow(2,INT) bits [33]
#   -H INT       use INT hash functions for Bloom filter [4]
#   -d FILE      dump hash table to FILE [null]
#   -E           skip error correction
#   -R           refine bfc-corrected reads
#   -r FILE      restore hash table from FILE [null]
#   -w INT       no more than 5 ec or 2 highQ ec in INT-bp window [10]
#   -c INT       min k-mer coverage [3]
#   -Q           force FASTA output
#   -1           drop reads containing unique k-mers
#   -v           show version number
#   -h           show command line help
