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
Reference Directory = $PIPELINE_REF
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
proceduresDir=$PIPELINE_HOME/procedures
dataDir=$PIPELINE_HOME/$subset
paramDir=$PIPELINE_HOME/$subset/model/$experiment/param/$parameters

format_status "Running Error Model Script"

# 
# Data Transfer
# 

# A Really Long & Dirty Conditional
if  [ "$experiment" = "bayeshammer" ] || \
    [ "$experiment" = "bfc" ]         || \
    [ "$experiment" = "blessec" ]     || \
    [ "$experiment" = "bloocoo" ]     || \
    [ "$experiment" = "decgpu" ]      || \
    [ "$experiment" = "karect" ]      || \
    [ "$experiment" = "kgem" ]        || \
    [ "$experiment" = "lighter" ]     || \
    [ "$experiment" = "musket" ]      || \
    [ "$experiment" = "quorum" ]      || \
    [ "$experiment" = "rcorrector" ]  || \
    [ "$experiment" = "nomodel" ]; then 


    if [ "$experiment" = "nomodel" ]; then

        format_status "Copying Read FASTQ Files to Modeled Directory..."
        files=$(echo $(ls $dataDir/fastq/split/$fileprefix.$subset.$condition.*.fastq))

        for file in $files
            # In Parallel
            do (
                # Extract Read Group to Pass Through
                suffix=$(echo "$file" | sed "s|$dataDir/fastq/split/$fileprefix.$subset.$condition.||")
                readgroup=$(echo "$suffix" | sed "s|.fastq$||")
                substate="$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup:ERRORMODEL:1"
                
                # Run Command
                if !(has_state $substate); then

                    # No Model, Copy Data to Modeled Directory
                    # Define Command
                    call="cp $file $paramDir/modeled/$fileprefix.$subset.$condition.$experiment.$parameters.$readgroup.fastq"
                    # Print & Call
                    format_status "Command:\n$call"
                    eval $call

                    # Check for failed parallel call
                    put_state $? $substate

                fi
            ) &
        done
        wait # Prevent Premature Exiting of Script
        
    else

        # State Check - Run Block if it Has Not Already Been Executed Successfully
        state="$fileprefix.$subset.$condition.$experiment.$parameters:ERRORMODEL:1"
        if !(has_state $state); then

            # Copy Data to Pre-Align Directory For Error Models
            format_status "Copying Read FASTQ Files to Pre-Alignment Directory..."
            # Define Command
            call="cp $dataDir/fastq/split/$fileprefix.$subset.$condition.*.fastq $paramDir/pre-align/fastq/"
            # Print & Call
            format_status "Command:\n$call"
            eval $call
    
            # Update State on Exit
            put_state $? $state

        fi

    fi

    format_status "FASTQ Copy Complete"

fi

#
# Delegate Args & Call Experiment Script
#

# State Management Delegated to Substates
case "$experiment" in

    "bayeshammer")

        format_status "Running Bayes Hammer"
        # Retrieve Files to Process
        files=$(echo $(ls $paramDir/pre-align/fastq/$fileprefix.$subset.$condition.*.fastq))
        for file in $files
            do (
                # Extract Read Group to Pass Through
                suffix=$(echo "$file" | sed "s|$paramDir/pre-align/fastq/$fileprefix.$subset.$condition.||")
                readgroup=$(echo "$suffix" | sed "s|.fastq$||")
                # Run Error Model by Readgroup
                # Define Command
                call="source $proceduresDir/models/bayeshammer.sh -f=$fileprefix -s=$subset -c=$condition -g=$readgroup -x=$experiment -p=$parameters -n=$ncores -m=$memory"
                # Print & Call
                format_status "Command:\n$call"
                eval $call 
            )
        done
        format_status "Error Model Complete"
    ;;

    "bfc")

        format_status "Running BFC"
        # Retrieve Files to Process
        files=$(echo $(ls $paramDir/pre-align/fastq/$fileprefix.$subset.$condition.*.fastq))
        for file in $files
            do (
                # Extract Read Group to Pass Through
                suffix=$(echo "$file" | sed "s|$paramDir/pre-align/fastq/$fileprefix.$subset.$condition.||")
                readgroup=$(echo "$suffix" | sed "s|.fastq$||")
                # Run Error Model by Readgroup
                # Define Command
                call="source $proceduresDir/models/bfc.sh -f=$fileprefix -s=$subset -c=$condition -g=$readgroup -x=$experiment -p=$parameters -n=$ncores -m=$memory"
                # Print & Call
                format_status "Command:\n$call"
                eval $call 
            )
        done
        format_status "Error Model Complete"
    ;;

    "blessec")

        format_status "Running Bless-EC"
        # Retrieve Files to Process
        files=$(echo $(ls $paramDir/pre-align/fastq/$fileprefix.$subset.$condition.*.fastq))
        for file in $files
            do (
                # Extract Read Group to Pass Through
                suffix=$(echo "$file" | sed "s|$paramDir/pre-align/fastq/$fileprefix.$subset.$condition.||")
                readgroup=$(echo "$suffix" | sed "s|.fastq$||")
                # Run Error Model by Readgroup
                # Define Command
                call="source $proceduresDir/models/blessec.sh -f=$fileprefix -s=$subset -c=$condition -g=$readgroup -x=$experiment -p=$parameters -n=$ncores -m=$memory"
                # Print & Call
                format_status "Command:\n$call"
                eval $call 
            )
        done
        format_status "Error Model Complete"
    ;;

    "bloocoo")

        format_status "Running Bloocoo"
        # Retrieve Files to Process
        files=$(echo $(ls $paramDir/pre-align/fastq/$fileprefix.$subset.$condition.*.fastq))
        for file in $files
            do (
                # Extract Read Group to Pass Through
                suffix=$(echo "$file" | sed "s|$paramDir/pre-align/fastq/$fileprefix.$subset.$condition.||")
                readgroup=$(echo "$suffix" | sed "s|.fastq$||")
                # Run Error Model by Readgroup
                # Define Command
                call="source $proceduresDir/models/bloocoo.sh -f=$fileprefix -s=$subset -c=$condition -g=$readgroup -x=$experiment -p=$parameters -n=$ncores -m=$memory"
                # Print & Call
                format_status "Command:\n$call"
                eval $call 
            )
        done
        format_status "Error Model Complete"
    ;;
    "decgpu")

        format_status "Running Dec-GPU"
         # Retrieve Files to Process
        files=$(echo $(ls $paramDir/pre-align/fastq/$fileprefix.$subset.$condition.*.fastq))
        for file in $files
            do (
                # Extract Read Group to Pass Through
                suffix=$(echo "$file" | sed "s|$paramDir/pre-align/fastq/$fileprefix.$subset.$condition.||")
                readgroup=$(echo "$suffix" | sed "s|.fastq$||")
                # Run Error Model by Readgroup
                # Define Command
                call="source $proceduresDir/models/decgpu.sh -f=$fileprefix -s=$subset -c=$condition -g=$readgroup -x=$experiment -p=$parameters -n=$ncores -m=$memory"
                # Print & Call
                format_status "Command:\n$call"
                eval $call 
            )
        done
        format_status "Error Model Complete"
    ;;

    "karect")

        format_status "Running Karect"
        # Retrieve Files to Process
        files=$(echo $(ls $paramDir/pre-align/fastq/$fileprefix.$subset.$condition.*.fastq))
        for file in $files
            do (
                # Extract Read Group to Pass Through
                suffix=$(echo "$file" | sed "s|$paramDir/pre-align/fastq/$fileprefix.$subset.$condition.||")
                readgroup=$(echo "$suffix" | sed "s|.fastq$||")
                # Run Error Model by Readgroup
                # Define Command
                call="source $proceduresDir/models/karect.sh -f=$fileprefix -s=$subset -c=$condition -g=$readgroup -x=$experiment -p=$parameters -n=$ncores -m=$memory"
                # Print & Call
                format_status "Command:\n$call"
                eval $call 
            )
        done
        format_status "Error Model Complete"
    ;;

    "kgem")

        format_status "Running KGEM"
        # Retrieve Files to Process
        files=$(echo $(ls $paramDir/pre-align/fastq/$fileprefix.$subset.$condition.*.fastq))
        for file in $files
            do (
                # Extract Read Group to Pass Through
                suffix=$(echo "$file" | sed "s|$paramDir/pre-align/fastq/$fileprefix.$subset.$condition.||")
                readgroup=$(echo "$suffix" | sed "s|.fastq$||")
                # Run Error Model by Readgroup
                # java Xmx$memory -jar $ERIF
                # java Xmx$memory -jar $KGEM
                # Define Command
                call="source $proceduresDir/models/kgem.sh -f=$fileprefix -s=$subset -c=$condition -g=$readgroup -x=$experiment -p=$parameters -n=$ncores -m=$memory"
                # Print & Call
                format_status "Command:\n$call"
                eval $call 
            )
        done
        format_status "Error Model Complete"
    ;;

    "lighter")

        format_status "Running Lighter"
        # Retrieve Files to Process
        files=$(echo $(ls $paramDir/pre-align/fastq/$fileprefix.$subset.$condition.*.fastq))
        for file in $files
            do (
                # Extract Read Group to Pass Through
                suffix=$(echo "$file" | sed "s|$paramDir/pre-align/fastq/$fileprefix.$subset.$condition.||")
                readgroup=$(echo "$suffix" | sed "s|.fastq$||")
                # Run Error Model by Readgroup
                # Define Command
                call="source $proceduresDir/models/lighter.sh -f=$fileprefix -s=$subset -c=$condition -g=$readgroup -x=$experiment -p=$parameters -n=$ncores -m=$memory"
                # Print & Call
                format_status "Command:\n$call"
                eval $call 
            )
        done
        format_status "Error Model Complete"
    ;;

    "musket")

        format_status "Running Musket"
         # Retrieve Files to Process
        files=$(echo $(ls $paramDir/pre-align/fastq/$fileprefix.$subset.$condition.*.fastq))
        for file in $files
            do (
                # Extract Read Group to Pass Through
                suffix=$(echo "$file" | sed "s|$paramDir/pre-align/fastq/$fileprefix.$subset.$condition.||")
                readgroup=$(echo "$suffix" | sed "s|.fastq$||")
                # Define Command
                call="source $proceduresDir/models/musket.sh -f=$fileprefix -s=$subset -c=$condition -g=$readgroup -x=$experiment -p=$parameters -n=$ncores -m=$memory"
                # Print & Call
                format_status "Command:\n$call"
                eval $call 
            )
        done
        format_status "Error Model Complete"
    ;;

    "quorum")
        
        format_status "Running Quorum"
         # Retrieve Files to Process
        files=$(echo $(ls $paramDir/pre-align/fastq/$fileprefix.$subset.$condition.*.fastq))
        for file in $files
            do (
                # Extract Read Group to Pass Through
                suffix=$(echo "$file" | sed "s|$paramDir/pre-align/fastq/$fileprefix.$subset.$condition.||")
                readgroup=$(echo "$suffix" | sed "s|.fastq$||")
                # Define Command
                call="source $proceduresDir/models/quorum.sh -f=$fileprefix -s=$subset -c=$condition -g=$readgroup -x=$experiment -p=$parameters -n=$ncores -m=$memory"
                # Print & Call
                format_status "Command:\n$call"
                eval $call 
            )
        done
        format_status "Error Model Complete"
    ;;

    "rcorrector")

        format_status "Running Rcorrector"
        # Retrieve Files to Process
        files=$(echo $(ls $paramDir/pre-align/fastq/$fileprefix.$subset.$condition.*.fastq))
        for file in $files
            do (
                # Extract Read Group to Pass Through
                suffix=$(echo "$file" | sed "s|$paramDir/pre-align/fastq/$fileprefix.$subset.$condition.||")
                readgroup=$(echo "$suffix" | sed "s|.fastq$||")
                # Run Error Model by Readgroup
                # perl $RCORRECTOR
                # Define Command
                call="source $proceduresDir/models/rcorrector.sh -f=$fileprefix -s=$subset -c=$condition -g=$readgroup -x=$experiment -p=$parameters -n=$ncores -m=$memory"
                # Print & Call
                format_status "Command:\n$call"
                eval $call 
            )
        done
        format_status "Error Model Complete"
    ;;

    "seecer")

        format_status "Running Seecer"
        # Retrieve Files to Process
        files=$(echo $(ls $paramDir/pre-align/fastq/$fileprefix.$subset.$condition.*.fastq))
        for file in $files
            do (
                # Extract Read Group to Pass Through
                suffix=$(echo "$file" | sed "s|$paramDir/pre-align/fastq/$fileprefix.$subset.$condition.||")
                readgroup=$(echo "$suffix" | sed "s|.fastq$||")
                # Run Error Model by Readgroup
                # source $SEECER
                # Define Command
                call="source $proceduresDir/models/seecer.sh -f=$fileprefix -s=$subset -c=$condition -g=$readgroup -x=$experiment -p=$parameters -n=$ncores -m=$memory"
                # Print & Call
                format_status "Command:\n$call"
                eval $call 
            )
        done
        format_status "Error Model Complete"
    ;;

    "shorah")

        format_status "Running SHoRAH"
        # Retrieve Files to Process
        files=$(echo $(ls $paramDir/pre-align/fastq/$fileprefix.$subset.$condition.*.fastq))
        for file in $files
            do (
                # Extract Read Group to Pass Through
                suffix=$(echo "$file" | sed "s|$paramDir/pre-align/fastq/$fileprefix.$subset.$condition.||")
                readgroup=$(echo "$suffix" | sed "s|.fastq$||")
                # Run Error Model by Readgroup
                # python $SHORAH
                # Define Command
                call="source $proceduresDir/models/shorah.sh -f=$fileprefix -s=$subset -c=$condition -g=$readgroup -x=$experiment -p=$parameters -n=$ncores -m=$memory"
                # Print & Call
                format_status "Command:\n$call"
                eval $call 
            )
        done
        format_status "Error Model Complete"
    ;;

    "nomodel")
        format_status "No Model Selected. Exiting."
    ;;

    "norealign")
        format_status "No Realignment Selected"
        format_status "Nothing to do. Exiting."
    ;;

    # Catch Any Invalid Error Models & Output Error
    *)

    format_status "Invalid Experiment (Error Model) Parameter: $experiment"
    ;;

esac

# Note State Completion is Managed by Fulfillment of Model Child States
format_status "Done"
