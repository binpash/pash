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
maxMemory=$((allocMemory * ncores))$allocSize

# Set Directories
dataDir=$PIPELINE_HOME/$subset
paramDir=$PIPELINE_HOME/$subset/model/$experiment/param/$parameters

# Test for Paired Ends
head -n 1000 $dataDir/fastq/split/$fileprefix.$subset.$condition.$readgroup.fastq | grep -qE "^@.*/1"
end1=$?
head -n 1000 $dataDir/fastq/split/$fileprefix.$subset.$condition.$readgroup.fastq | grep -qE "^@.*/2"
end2=$?

# Is Interleaved?
paired=""

#
# Paired End Input
#

if [ $end1 ] && [ $end2 ]; then

    #
    # FASTQ File Pre-Processing
    #

    # State Check - Run Block if it Has Not Already Been Executed Successfully
    state="$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup:QUORUM:1"
    if !(has_state $state); then

        format_status "Interleaved Paired-End Detected (/1 = $end1, /2 = $end2)"
        format_status "Splitting Paired End FASTQ"
        # Define Command
        call="$PYTHON $PIPELINE_HOME/utils/paired-end-to-single-ends.py \
        -i $dataDir/fastq/split/$fileprefix.$subset.$condition.$readgroup.fastq \
        -o $dataDir/fastq/split/unpaired/$fileprefix.$subset.$condition.$readgroup"
        # Print & Call
        format_status "Command:\n$call"
        eval $call
        # Update State on Exit
        status=$?
        put_state $status $state

    fi

    # 
    # Run Quorum - Paired End
    # 

    # State Check - Run Block if it Has Not Already Been Executed Successfully
    state="$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup:QUORUM:2"
    if !(has_state $state); then

        # Gather Unpaired Reads
        input=$dataDir/fastq/split/unpaired/$fileprefix.$subset.$condition.$readgroup.*.fastq

        if [ "$parameters" = "default" ]; then

            #
            # Default Parameters
            #
            
            format_status "Running Quorum - $parameters Parameters"
            mkdir -p $paramDir/modeled/unpaired
            # Define Command
            call="$QUORUM \
            $input \
            --prefix $paramDir/modeled/unpaired/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup \
            -t $ncores \
            --size 48000000000 \
            --no-discard \
            --min-q-char 33 \
            --paired-files \
            --debug"
            # Print & Call
            format_status "Command:\n$call"
            eval $call

        elif [ "$parameters" = "custom" ]; then

            #
            # Custom Parameters
            #
            
            format_status "Running Quorum - $parameters Parameters"
            mkdir -p $paramDir/modeled/unpaired
            # Define Command
            call="$QUORUM \
            $input \
            --prefix $paramDir/modeled/unpaired/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup \
            -t $ncores \
            --size 48000000000 \
            --no-discard \
            --min-q-char 33 \
            --paired-files \
            --debug"
            # Print & Call
            format_status "Command:\n$call"
            eval $call

        fi

        # Update State on Exit
        status=$?
        put_state $status $state
        format_status "Quorum ($parameters $readgroup) Complete"

    fi

    # 
    # Interleave Split Pairs
    # 

    # State Check - Run Block if it Has Not Already Been Executed Successfully
    state="$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup:QUORUM:3"
    if !(has_state $state); then

        format_status "Merging Split FASTQ Pairs"
        # Define Command
        call="$PYTHON $PIPELINE_HOME/utils/single-ends-to-paired-end.py \
        -1 $paramDir/modeled/unpaired/$fileprefix.$subset.$condition.$readgroup.1.fastq \
        -2 $paramDir/modeled/unpaired/$fileprefix.$subset.$condition.$readgroup.2.fastq \
        -o $paramDir/modeled/$fileprefix.$subset.$condition.$readgroup"
        # Print & Call
        format_status "Command:\n$call"
        eval $call
        # Update State on Exit
        status=$?
        put_state $status $state

        return $status

    fi

else

    # 
    # Run Quorum - Single End
    # 

    # State Check - Run Block if it Has Not Already Been Executed Successfully
    state="$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup:QUORUM:1"
    if !(has_state $state); then

        format_status "Single-End Detected  (/1 = $end1, /2 = $end2)"
        input=$dataDir/fastq/split/$fileprefix.$subset.$condition.$readgroup.fastq

        if [ "$parameters" = "default" ]; then

            #
            # Default Parameters
            #
            
            format_status "Running Quorum - $parameters Parameters"
            # Define Command
            call="$QUORUM \
            $input \
            --prefix $paramDir/modeled/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup \
            -t $ncores \
            --size 43000000000 \
            --no-discard \
            --min-q-char 33 \
            --debug"
            # Print & Call
            format_status "Command:\n$call"
            eval $call

        elif [ "$parameters" = "custom" ]; then

            #
            # Custom Parameters
            #
            
            format_status "Running Quorum - $parameters Parameters"
            # Define Command
            call="$QUORUM \
            $input \
            --prefix $paramDir/modeled/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup \
            -t $ncores \
            --size 43000000000 \
            --no-discard \
            --min-q-char 33 \
            --debug"
            # Print & Call
            format_status "Command:\n$call"
            eval $call

        fi

        # Update State on Exit
        status=$?
        put_state $status $state
        format_status "Quorum ($parameters $readgroup) Complete"

        return $status

    fi

fi

#
# Quorum Arguments
#

# -s, --size              Mer database size (default 200M)
# -t, --threads           Number of threads (default number of cpus)
# -p, --prefix            Output prefix (default quorum_corrected)
# -k, --kmer-len          Kmer length (default 24)
# -q, --min-q-char        Minimum quality char. Usually 33 or 64 (autodetect)
# -m, --min-quality       Minimum above -q for high quality base (5)
# -w, --window            Window size for trimming
# -e, --error             Maximum number of errors in a window
#     --min-count         Minimum count for a k-mer to be good
#     --skip              Number of bases to skip to find anchor kmer
#     --anchor            Numer of good kmer in a row for anchor
#     --anchor-count      Minimum count for an anchor kmer
#     --contaminant       Contaminant sequences
#     --trim-contaminant  Trim sequences with contaminant mers
# -d, --no-discard        Do not discard reads, output a single N (false)
# -P, --paired-files      Preserve mate pairs in two files
#     --homo-trim         Trim homo-polymer on 3' end
#     --debug             Display debugging information
#     --version           Display version
# -h, --help              This message
