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

        # n Reads Per GB or Memory
        -r=*|--reads=*)
        readsOpt="${i#*=}"
        shift
        ;;

        # Read Group Regex
    	-g=*|--regex=*)
    	regexOpt="${i#*=}"
    	shift
    	;;

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

        # Trigger Debugging Available in Tools
        --debug)
        debugOpt=true
        ;;

        # Pass Through, Not Applying Mark Duplicates
        --pass)
        passOpt=true
        ;;

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
readsDef=150000
regexDef=""
ncoresDef="16"
memoryDef="6G"
passDef=false
debugDef=false
cleanDef=false

# Set Optional Values
reads=${readsOpt:-$readsDef}
regex=${regexOpt:-$regexDef}
ncores=${ncoresOpt:-$ncoresDef}
memory=${memoryOpt:-$memoryDef}
pass=${passOpt:-$passDef}
debug=${debugOpt:-$debugDef}
clean=${cleanOpt:-$cleanDef}

# Get Max Allowable Memory
allocMemory=${memory//[GgMmKk]/}
allocSize=${memory//[0-9]/}
maxMemory=$((allocMemory * ncores))$allocSize

# Max Reads in RAM - 200,000 per GB
maxReads=$((allocMemory * $reads))

printf "\nPARAMETERS: 
Picard Directory    = $PICARD
Data File Prefix    = $fileprefix
Data Subset         = $subset
Condition           = $condition
Experiment          = $experiment
Parameter Set       = $parameters
Read Group Regex    = $regex
Memory              = $memory
Cores               = $ncores
Max Memory          = $maxMemory
Max Reads in Memory = $maxReads
Debug               = $debug
Pass                = $pass
Do Cleanup          = $clean
\n"

# Set Directories
dataDir=$PIPELINE_HOME/$subset
paramDir=$PIPELINE_HOME/$subset/model/$experiment/param/$parameters
tmpDir=$PIPELINE_HOME/$subset/tmp

# Set Readgroup Regex
if [ "$regex" = "" ]; then
    rgRegex=""
else
    rgRegex="READ_NAME_REGEX=\"$regex\""
fi

# Tool Specific Debugging - Picard
verbosity="INFO"

if $debug; then 
    verbosity="DEBUG"
fi

format_status "Running Mark Duplicates Script"

# If Norealignment, Get Data From Download Directory & Skip Sorting
if [ "$experiment" = "norealign" ] && [ $pass == false ]; then

    #
    # Mark Duplicates
    #

    # State Check - Run Block if it Has Not Already Been Executed Successfully
    state="$fileprefix.$subset.$condition.$experiment.$parameters:MARKDUPLICATES:1"
    if !(has_state $state); then

        format_status "MarkDuplicates Start"

        # Define Command
        call="java -Xmx$maxMemory \
        -Djava.io.tmpdir=$tmpDir \
        -jar $PICARD MarkDuplicates \
        I=$dataDir/downloaded/$fileprefix.$subset.$condition.bam \
        O=$paramDir/markdup/$fileprefix.$subset.$condition.$experiment.$parameters.bam \
        M=$paramDir/markdup/log_marked_duplicates_metrics_$condition.txt \
        MAX_RECORDS_IN_RAM=$maxReads \
        PG=null \
        TMP_DIR=$tmpDir \
        VERBOSITY=$verbosity"

        # Print & Call
        format_status "Command:\n$call $rgRegex"
        eval $call $rgRegex

        # Update State on Exit
        put_state $? $state
        format_status "Mark Duplicates Complete"

    fi

    #
    # Create New BAM Index
    #

    # State Check - Run Block if it Has Not Already Been Executed Successfully
    state="$fileprefix.$subset.$condition.$experiment.$parameters:MARKDUPLICATES:2"
    if !(has_state $state); then

        format_status "Indexing BAM"
        # Define Command
        call="samtools index $paramDir/markdup/$fileprefix.$subset.$condition.$experiment.$parameters.bam"
        # Print & Call
        format_status "Command:\n$call"
        eval $call

        # Update State on Exit
        put_state $? $state
        format_status "BAM Indexing Complete"

    fi

# Otherwise, Grab Data From Merged Alignment Directory & Sort
else

    # 
    # Sort BAM
    # 

    # State Check - Run Block if it Has Not Already Been Executed Successfully
    state="$fileprefix.$subset.$condition.$experiment.$parameters:MARKDUPLICATES:1"
    if !(has_state $state); then

        format_status "Sorting BAM"
        # Define Command
        call="samtools sort -m $memory -@ $ncores -T $tmpDir $paramDir/merged/$fileprefix.$subset.$condition.$experiment.$parameters.bam -o $paramDir/merged/$fileprefix.$subset.$condition.$experiment.$parameters.sorted.bam"
        # Print & Call
        format_status "Command:\n$call"
        eval $call

        # Update State on Exit
        put_state $? $state
        format_status "Sort BAM Complete"

    fi

    #
    # Mark Duplicates
    #

    # If Passing, Copy File to Output Directory
    if $pass; then 

        # State Check - Run Block if it Has Not Already Been Executed Successfully
        state="$fileprefix.$subset.$condition.$experiment.$parameters:MARKDUPLICATES:2"
        if !(has_state $state); then

            format_status "Passing MarkDuplicates"
            # Define Command
            call="cp $paramDir/merged/$fileprefix.$subset.$condition.$experiment.$parameters.sorted.bam $paramDir/markdup/$fileprefix.$subset.$condition.$experiment.$parameters.bam"
            # Print & Call
            format_status "Command:\n$call"
            eval $call

            put_state $? $state
            format_status "Mark Duplicates Pass Complete"

        fi

        #
        # Create New BAM Index
        #

        # State Check - Run Block if it Has Not Already Been Executed Successfully
        state="$fileprefix.$subset.$condition.$experiment.$parameters:MARKDUPLICATES:3"
        if !(has_state $state); then

            format_status "Indexing BAM Output"
            # Define Command
            call="samtools index $paramDir/markdup/$fileprefix.$subset.$condition.$experiment.$parameters.bam"
            # Print & Call
            format_status "Command:\n$call"
            eval $call

            # Update State on Exit
            put_state $? $state
            format_status "BAM Indexing Complete"

        fi

    # Run Mark Duplicates
    else

        # State Check - Run Block if it Has Not Already Been Executed Successfully
        state="$fileprefix.$subset.$condition.$experiment.$parameters:MARKDUPLICATES:3"
        if !(has_state $state); then

            format_status "MarkDuplicates Start"
            # Define Command
            call="java -Xmx$memory \
            -Djava.io.tmpdir=$tmpDir \
            -jar $PICARD MarkDuplicates \
            I=$paramDir/merged/$fileprefix.$subset.$condition.$experiment.$parameters.sorted.bam \
            O=$paramDir/markdup/$fileprefix.$subset.$condition.$experiment.$parameters.bam \
            M=$paramDir/markdup/log_marked_duplicates_metrics_$condition.txt \
            MAX_RECORDS_IN_RAM=$maxReads \
            PG=null \
            TMP_DIR=$tmpDir \
            VERBOSITY=$verbosity"
            # Print & Call
            format_status "Command:\n$call $rgRegex"
            eval $call $rgRegex

            # Update State on Exit
            put_state $? $state
            format_status "Mark Duplicates Complete"

        fi

        #
        # Create New BAM Index
        #

        # State Check - Run Block if it Has Not Already Been Executed Successfully
        state="$fileprefix.$subset.$condition.$experiment.$parameters:MARKDUPLICATES:4"
        if !(has_state $state); then

            format_status "Indexing BAM Output"
            # Define Command
            call="samtools index $paramDir/markdup/$fileprefix.$subset.$condition.$experiment.$parameters.bam"
            # Print & Call
            format_status "Command:\n$call"
            eval $call

            # Update State on Exit
            put_state $? $state
            format_status "BAM Indexing Complete"

        fi

    fi

fi

format_status "Done"

