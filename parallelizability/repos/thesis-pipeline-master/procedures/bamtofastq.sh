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
   
    # Optional Arguments With Defaults

        -n=*|--ncores=*)
        ncoresOpt="${i#*=}"
        shift # Number of Cores to Use
        ;;

        -m=*|--memory=*)
        memoryOpt="${i#*=}"
        shift # Per Core Memory Requirement
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
Memory              = $memory
Cores               = $ncores
Max Memory          = $maxMemory
Do Cleanup          = $clean
\n"

# Set Directories
dataDir=$PIPELINE_HOME/$subset
tmpDir=$PIPELINE_HOME/$subset/tmp

format_status "Running BAM to FastQ Script"

#
# Split BAM
#

# State Check - Run Block if it Has Not Already Been Executed Successfully
state="$fileprefix.$subset.$condition:BAMTOFASTQ:1"
if !(has_state $state); then

    format_status "Splitting Merged BAM"
    # Define Command
    call="samtools split $dataDir/downloaded/$fileprefix.$subset.$condition.bam -f $dataDir/downloaded/split/$fileprefix.$subset.$condition.%!.%."
    # Print & Call
    format_status "Command:\n$call"
    eval $call

    # Update State on Exit
    put_state $? $state

fi

#
# Sort BAMs
#

# State Management Delegated to Substates
format_status "Running Samtools Sort BAM"
# Retrieve Files
files=$(echo $(ls $dataDir/downloaded/split/$fileprefix.$subset.$condition.*.bam))

for file in $files
    # In Parallel
    do ( 
        # Get Read Group to Process
        suffix=$(echo "$file" | sed "s|$dataDir/downloaded/split/$fileprefix.$subset.$condition.||")
        readgroup=$(echo "$suffix" | sed "s|.bam$||")
        substate="$fileprefix.$subset.$condition.$readgroup:BAMTOFASTQ:2"
        
        # State Check - Run Block if it Has Not Already Been Executed Successfully
        if !(has_state $substate); then
            
            # Call Bam to FastQ
            # Define Command
            call="samtools sort -n -m $memory -@ $ncores $file > $dataDir/downloaded/split/sorted/$fileprefix.$subset.$condition.$readgroup.bam"
            # Print & Call
            format_status "Command:\n$call"
            eval $call

            # Check for failed parallel call
            put_state $? $substate

        fi

    )

done

#
# Bam to FastQ
#

# State Management Delegated to Substates
format_status "Running Samtools BAM to FastQ"
# Retrieve Files
files=$(echo $(ls $dataDir/downloaded/split/sorted/$fileprefix.$subset.$condition.*.bam))

for file in $files
    # In Parallel
    do ( 
        # Get Read Group to Process
        suffix=$(echo "$file" | sed "s|$dataDir/downloaded/split/sorted/$fileprefix.$subset.$condition.||")
        readgroup=$(echo "$suffix" | sed "s|.bam$||")
        substate="$fileprefix.$subset.$condition.$readgroup:BAMTOFASTQ:3"
        
        # State Check - Run Block if it Has Not Already Been Executed Successfully
        if !(has_state $substate); then
            
            # Call Bam to FastQ
            # Define Command
            call="samtools fastq -t -O $file > $dataDir/fastq/split/$fileprefix.$subset.$condition.$readgroup.fastq"
            # Print & Call
            format_status "Command:\n$call"
            eval $call

            # Check for failed parallel call
            put_state $? $substate

        fi

    )

done

# Run Cleanup
if $clean; then
    rm $dataDir/downloaded/split/$fileprefix.$subset.$condition.*.bam
    rm $dataDir/downloaded/split/sorted/$fileprefix.$subset.$condition.*.bam
fi

format_status "Done"

