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

    # Optional Flags

        # Directory Cleanup (Voids All Other Parameters)
        --clean)
        cleanOpt=true
        ;;

    # Invalid Argument Handler

        *)
        # invalid option
        printf "Invalid/Unused Parameter: $i"
        ;;
        
    esac
done

# Defaults if No Arguments Passed
ncoresDef="16"
memoryDef="6G"
cleanDef=false

# Set Optional Values
ncores=${ncoresOpt:-$ncoresDef}
memory=${memoryOpt:-$memoryDef}
clean=${cleanOpt:-$cleanDef}

# Get Max Allowable Memory
allocMemory=${memory//[GgMmKk]/}
allocSize=${memory//[0-9]/}
maxMemory=$((allocMemory * ncores))$allocSize
 
printf "\nPARAMETERS: 
Data File Prefix    = $fileprefix
Data Subset         = $subset
Condition           = $condition
Experiment          = $experiment
Parameter Set       = $parameters
Memory              = $memory
Cores               = $ncores
Max Memory          = $maxMemory
Do Cleanup          = $clean
\n"

# Set Directories
dataDir=$PIPELINE_HOME/$subset
paramDir=$PIPELINE_HOME/$subset/model/$experiment/param/$parameters

format_status "Running Merge BAMs Script"

# 
# Inserting Read Groups
# 

# State Management Delegated to Substates
format_status "Samtools AddReplaceRG"
# Retrieve Files to Process
files=$(echo $(ls $dataDir/downloaded/split/$fileprefix.$subset.$condition.*.bam))

for file in $files
    do (
        # Get Read Group to Process
        suffix=$(echo "$file" | sed "s|$dataDir/downloaded/split/$fileprefix.$subset.$condition.||")
        readgroup=$(echo "$suffix" | sed "s|.bam$||")
        substate="$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup:MERGEALIGNMENT:1"
        
        # Run Command
        if !(has_state $substate); then
            # Get Read Group Arguments to Pass to Samtools
            rgArgs=$(samtools view -H $file | grep '@RG' | awk -F '\t' '{for (i=2; i<=NF; i++) print $i}' | sed "s|[A-Z][A-Z]:[a-zA-Z0-9\.\-\:]*|-r &|g")
            # Check for failed parallel call
            put_state $? $substate
        fi

        substate="$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup:MERGEALIGNMENT:2"
        
        # Run Command
        if !(has_state $substate); then

            # Insert Read Groups into New BAM
            # Define Command
            call="samtools addreplacerg ${rgArgs[@]} $paramDir/post-align/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.sam > $paramDir/post-align/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.bam"
            # Print & Call
            format_status "Command:\n$call"
            eval $call

            # Check for failed parallel call
            put_state $? $substate

        fi
    )
done

# Notify
format_status "Samtools AddReplaceRG Complete"

# 
# Merge BAMs to Single BAM
# 

# State Check - Run Block if it Has Not Already Been Executed Successfully
state="$fileprefix.$subset.$condition.$experiment.$parameters:MERGEALIGNMENT:2"
if !(has_state $state); then

    # Retrieve Files to Process
    files=$(echo $(ls $paramDir/post-align/$fileprefix.$subset.$condition.$experiment.$parameters.*.bam))

    format_status "Samtools Merge"
    # Define Command
    call="samtools merge -r -f $paramDir/merged/$fileprefix.$subset.$condition.$experiment.$parameters.bam $files"
    # Print & Call
    format_status "Command:\n$call"
    eval $call

    # Update State on Exit
    put_state $? $state
    format_status "Samtools Merge Complete"

fi

format_status "Done"
