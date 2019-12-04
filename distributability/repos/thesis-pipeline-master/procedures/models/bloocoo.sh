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

# Workaround to Bloocoo's Funky Output Scheme (See State BLOOCOO:2)
corrected="_corrected"

# Set Directories
proceduresDir=$PIPELINE_HOME/procedures
dataDir=$PIPELINE_HOME/$subset
paramDir=$PIPELINE_HOME/$subset/model/$experiment/param/$parameters

#
# Convert FastQ to Fasta & Qual
#

# State Check - Run Block if it Has Not Already Been Executed Successfully
state="$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup:BLOOCOO:1"
if !(has_state $state); then

        format_status "Splitting FastQ to Fasta & Qual"
        # Define Command
        call="$PYTHON $PIPELINE_HOME/utils/fastq-to-fasta-qual.py \
        -i $dataDir/fastq/split/$fileprefix.$subset.$condition.$readgroup.fastq \
        -o $dataDir/fastq/split/unmerged/$fileprefix.$subset.$condition.$readgroup"
        # Print & Call
        format_status "Command:\n$call"
        eval $call

        # Update State on Exit
        status=$?
        put_state $status $state

fi

# 
# Run Bloocoo
# 

# State Check - Run Block if it Has Not Already Been Executed Successfully
state="$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup:BLOOCOO:2"
if !(has_state $state); then

    if [ "$parameters" = "default" ]; then
        
        #
        # Default Parameters
        #

        format_status "Running Bloocoo - $parameters Parameters"
        # Define Command
        call="$BLOOCOO \
        -file $dataDir/fastq/split/unmerged/$fileprefix.$subset.$condition.$readgroup.fasta \
        -nb-cores $ncores"
        # Print & Call
        format_status "Command:\n$call"
        eval $call

    elif [ "$parameters" = "custom" ]; then
        
        format_status "Running Bloocoo - Custom Parameters"
        # Define Command
        call="$BLOOCOO \
        -file $dataDir/fastq/split/unmerged/$fileprefix.$subset.$condition.$readgroup.fasta \
        -nb-cores $ncores \
        -slow \
        -high-precision"
        # Print & Call
        format_status "Command:\n$call"
        eval $call
        

    fi

    # Update State on Exit
    status=$?
    put_state $status $state

fi

#
# Merge Fasta & Qual Back to FastQ & Output to Model Directory
#

# State Check - Run Block if it Has Not Already Been Executed Successfully
state="$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup:BLOOCOO:3"
if !(has_state $state); then

        format_status "Merging Fasta & Qual to FastQ"
        # Define Command
        call="$PYTHON $PIPELINE_HOME/utils/fasta-qual-to-fastq.py \
        -f $proceduresDir/$fileprefix.$subset.$condition.$readgroup$corrected.fasta \
        -q $dataDir/fastq/split/unmerged/$fileprefix.$subset.$condition.$readgroup.qual \
        -o $paramDir/modeled/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup"
        # Print & Call
        format_status "Command:\n$call"
        eval $call

        # Update State on Exit
        status=$?
        put_state $status $state
        
fi

#
# Cleanup
#

# State Check - Run Block if it Has Not Already Been Executed Successfully
state="$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup:BLOOCOO:4"
if !(has_state $state); then

        format_status "Cleaning Up Fasta & Qual Files"
        # Define Command
        call="rm $proceduresDir/$fileprefix.$subset.$condition.$readgroup$corrected.fasta"
        # Print & Call
        format_status "Command:\n$call"
        eval $call

        # Update State on Exit
        status=$?
        put_state $status $state

fi

#
# Bloocoo Arguments
#

# -file : the read file name, can be fasta, fastq, gzipped or not.
# -kmer-size : the k-mer size (typically ~31)
# -abundance-min : the minimal abundance threshold defining solid k-mers (typically  between 3 and 6, but depends on the read depth, you can also use 'auto' and it is automatically inferred from the data)
# -nb-cores  : number of threads used
# -high-recall  :  correct more errors but can also introduce more mistakes
# -slow : slower modes with more pass, but better correction
# -high-precision :  correct safely, correct less errors but introduce less mistakes
# -ion : (experimental) mode for correcting indels present in  ion torrent reads