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
# Unpair FASTQ File
#

# State Check - Run Block if it Has Not Already Been Executed Successfully
state="$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup:LIGHTER:1"
if !(has_state $state); then

     format_status "Splitting Paired End FASTQ to Single End"
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
# Run Lighter
# 

# State Check - Run Block if it Has Not Already Been Executed Successfully
state="$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup:LIGHTER:2"
if !(has_state $state); then

    if [ "$parameters" = "default" ]; then
        
        #
        # Default Parameters
        #

        format_status "Running Lighter - Default Parameters"
        # Define Command
        call="lighter \
        -od $paramDir/modeled \
        -r $dataDir/fastq/split/unpaired/$fileprefix.$subset.$condition.$readgroup.1.fastq \
        -r $dataDir/fastq/split/unpaired/$fileprefix.$subset.$condition.$readgroup.2.fastq \
        -t $ncores \
        -k 23 3250000000 0.225"
        # Print & Call
        format_status "Command:\n$call"
        eval $call

    elif [ "$parameters" = "custom" ]; then

        #
        # Custom Parameters
        #

        format_status "Running Lighter - Custom Parameters"
        # Define Command
        call="lighter \
        -od $paramDir/modeled \
        -r $dataDir/fastq/split/unpaired/$fileprefix.$subset.$condition.$readgroup.1.fastq \
        -r $dataDir/fastq/split/unpaired/$fileprefix.$subset.$condition.$readgroup.2.fastq \
        -t $ncores \
        -k 23 3250000000 0.1"
        # Print & Call
        format_status "Command:\n$call"
        eval $call

    fi
    
    # Update State on Exit
    status=$?
    put_state $status $state

fi

# 
# Merge Corrected Reads
# 

# State Check - Run Block if it Has Not Already Been Executed Successfully
state="$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup:LIGHTER:3"
if !(has_state $state); then

    format_status "Merging Single End FASTQ to Interleaved Paired End"
    # Define Command
    call="$PYTHON $PIPELINE_HOME/utils/single-ends-to-paired-end.py \
    -1 $paramDir/modeled/$fileprefix.$subset.$condition.$readgroup.1.cor.fq \
    -2 $paramDir/modeled/$fileprefix.$subset.$condition.$readgroup.2.cor.fq \
    -o $paramDir/modeled/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup"
    # Print & Call
    format_status "Command:\n$call"
    eval $call

    # Update State on Exit
    status=$?
    put_state $status $state
    format_status "Lighter ($parameters $readgroup) Complete"
    return $status

fi

# lighter
# Usage: ./lighter [OPTIONS]
# OPTIONS:
# Required parameters:
# 	-r seq_file: seq_file is the path to the sequence file. Can use multiple -r to specifiy multiple sequence files
# 	             The file can be fasta and fastq, and can be gzip'ed with extension *.gz.
# 	             When the input file is *.gz, the corresponding output file will also be gzip'ed.
# 	-k kmer_length genome_size alpha: (see README for information on setting alpha)
# 					or
# 	-K kmer_length genom_size: in this case, the genome size should be relative accurate.
# Other parameters:
# 	-od output_file_directory: (default: ./)
# 	-t num_of_threads: number of threads to use (default: 1)
# 	-maxcor INT: the maximum number of corrections within a 20bp window (default: 4)
# 	-trim: allow trimming (default: false)
# 	-discard: discard unfixable reads. Will LOSE paired-end matching when discarding (default: false)
# 	-noQual: ignore the quality socre (default: false)
# 	-newQual ascii_quality_score: set the quality for the bases corrected to the specified score (default: not used)
# 	-zlib compress_level: set the compression level(0-9) of gzip (default: 1)
# 	-h: print the help message and quit
# 	-v: print the version information and quit