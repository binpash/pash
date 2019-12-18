#! /usr/bin/bash

#
# A Positively Disgusting Submit Script Generator
#

timestamp=$(date +"%m-%d-%y_%H꞉%M")

# Exit on First Error - to Prevent Invalid File Modifications
# set -o errexit

# Assign Arguments
for i in "$@"
    do case $i in

    # Standard Arguments
        -b=*|--binary=*)
        binary="${i#*=}"
        shift # Shell Script Module to Call
        ;;
        -f=*|--fileprefix=*)
        fileprefix="${i#*=}"
        shift # Access & Write Files With This Prefix
        ;;
        -s=*|--subset=*)
        subset="${i#*=}"
        shift # Access & Write Files With This Subset
        ;;
        -c=*|--condition=*)
        condition="${i#*=}"
        shift # Access & Write Files With This Condition
        ;;
        -x=*|--experiment=*)
        experiment="${i#*=}"
        shift # Access & Write Files With This Experiment
        ;;
        -p=*|--parameters=*)
        parameters="${i#*=}"
        shift # Access & Write Files With This Parameter Set
        ;;
        -q=*|--qualitymodel=*)
        qualitymodel="${i#*=}"
        shift # Access & Write Files With This Quality Model
        ;;

    # Optional Arguments With Defaults

        -r=*|--reads=*)
        readsOpt="${i#*=}"
        shift # n Reads Per GB or Memory
        ;;
        -n=*|--ncores=*)
        ncoresOpt="${i#*=}"
        shift # Number of Cores to Use
        ;;
        -m=*|--memory=*)
        memoryOpt="${i#*=}"
        shift # Per Core Memory Requirement
        ;;
        -d=*|--debug=*)
        debugOpt="${i#*=}"
        shift # Trigger Debugging Available in Tools
        ;;

    # Tool Specific Arguments

        -A=*|--align=*)
        alignOpt="${i#*=}"
        shift # BWA Alignment Type
        ;;
        -C=*|--contamination=*)
        contaminationOpt="${i#*=}"
        shift # BWA Alignment Type
        ;;

    # Condor Only Arguments

        --notify=*)
        notifyOpt="${i#*=}"
        shift # Force job duration
        ;;
        --autosubmit=*)
        autosubmitOpt="${i#*=}"
        shift # Force job duration
        ;;
        --force=*)
        forceOpt="${i#*=}"
        shift # Force job duration
        ;;
        --priority=*)
        priorityOpt="${i#*=}"
        shift # Condor Priority (-20 to 20)
        ;;
        --universe=*)
        universeOpt="${i#*=}"
        shift # Condor universe
        ;;
        --suspendable=*)
        suspendableOpt="${i#*=}"
        shift # Is job suspendable
        ;;
        --getenv=*)
        getenvOpt="${i#*=}"
        shift # Get Shell Environment
        ;;
        --maxtime=*)
        maxtimeOpt="${i#*=}"
        shift # Get Shell Environment
        ;;

    # Invalid Argument Handler

        *)
        # invalid option
        printf "Invalid/Unused Parameter: $i"
        ;;
        
    esac
done

#
# Setup
# 

# Defaults if No Arguments Passed
ncoresDef="12"
memoryDef="8G"
readsDef=150000
debugDef=false

# Set Optional Values
ncores=${ncoresOpt:-$ncoresDef}
memory=${memoryOpt:-$memoryDef}
reads=${readsOpt:-$readsDef}
debug=${debugOpt:-$debugDef}

# Condor Defaults
autosubmitDef=true
forceDef=false
priorityDef=0
universeDef="vanilla"
suspendableDef=false
getenvDef=true
notifyDef=""
maxtimeDef="129600"

# Condor Set Optional Values
autosubmit=${autosubmitOpt:-$autosubmitDef}
force=${forceOpt:-$forceDef}
priority=${priorityOpt:-$priorityDef}
universe=${universeOpt:-$universeDef}
suspendable=${suspendableOpt:-$suspendableDef}
getenv=${getenvOpt:-$getenvDef}
notify=${notifyOpt:-$notifyDef}
maxtime=${maxtimeOpt:-$maxtimeDef}

# Get Max Allowable Memory
allocMemory=${memory//[GgMmKk]/}
allocSize=${memory//[0-9]/}
module=${binary//.sh/}
maxMemory=$((allocMemory * ncores))$allocSize
maxMemoryMB=$((allocMemory * ncores * 1024))

# Generate Argument String & Namespace File Handles
subfile=$timestamp.$module.$fileprefix.$subset
logfile=$timestamp.$module.$fileprefix.$subset
errfile=$timestamp.$module.$fileprefix.$subset
args=""

#
# Namespacing & Argument Concatenation
#

# By Prefix
if [[ ! -z "$fileprefix" ]]; then
    args=$args"-f=$fileprefix "
fi

# By Condition
if [[ ! -z "$condition" ]]; then
    subfile=$subfile.$condition
    logfile=$logfile.$condition
    errfile=$errfile.$condition
    args=$args"-x=$condition "
fi

# By Experiment
if [[ ! -z "$experiment" ]]; then
    subfile=$subfile.$experiment
    logfile=$logfile.$experiment
    errfile=$errfile.$experiment
    args=$args"-x=$experiment "
fi

# By Parameters
if [[ ! -z "$parameters" ]]; then
    subfile=$subfile.$parameters
    logfile=$logfile.$parameters
    errfile=$errfile.$parameters
    args=$args"-p=$parameters "
fi

# By Quality Model
if [[ ! -z $qualitymodel ]]; then
    subfile=$subfile.$qualitymodel
    logfile=$logfile.$qualitymodel
    errfile=$errfile.$qualitymodel
    args=$args"-q=$qualitymodel "
fi

# Append Suffix
logfile=$logfile.log
errfile=$errfile.err
subfile=$subfile.sub

#
# Append Additional Args to Args String
#

# Memory
if [ $memory ]; then
    args=$args"-m=$memory "
fi

# n Cores
if [ $ncores ]; then
    args=$args"-n=$ncores "
fi

# Reads per GB
if [ $reads ]; then
    args=$args"-r=$reads "
fi

# Debug Script
if [ $debug ]; then
    args=$args"-d=$debug "
fi
# BWA Alignment Type ('mem' or 'bsws')
if [ $align ]; then
    args=$args+"-a=$align "
fi
# ContEst Contamination Proportion (e.g. 0.001)
if [ $contamination ]; then
    args=$args+"-C=$contamination "
fi

#
# Generate Submit Script
#

format_status "Compiling Submit Script"

# Submit Header
header="####################################\n"

# Submit Boiler
submit="$header\n\n# Job Details\nexecutable = procedures/$binary\narguments = $args\nuniverse = $universe\npriority = $priority\n\n# Resource Requirements\nrequest_cpus = $ncores\nrequest_memory = $maxMemory\nimage_size = $maxMemory\nrank = Memory >= $maxMemoryMB\n\n# Logging\nlog = $HOME/logs/condor_jobs.log\noutput = submits/logs/$logfile\nerror = submits/logs/$errfile\n\n# Additional Arguments\n+MaxExecutionTime = $maxtime\nshould_transfer_files = IF_NEEDED\nwhen_to_transfer_output = ON_EXIT\n\n# Compiled Optional Arguments\n"

# Append Getenv
if $getenv; then
    submit=$submit"getenv = True\n"
fi

# Append Notification
if [[ ! -z $notify ]]; then
    submit=$submit"notification = Complete\nnotify_user = $notify\n"
fi

# Append Force
if $force; then
    submit=$submit"concurrency_limits = WEEK_LONG_JOBS\n"
fi

# Allow Job Suspension
if $suspendable; then
    submit=$submit"+IsSuspensionJob = True\n"
fi

# Queue
submit=$submit"queue\n"
submit=$submit$header

#
# Write Submit
#

echo -e $submit > $PIPELINE_HOME/submits/$subfile
echo -e "Submit Script Written to: $PIPELINE_HOME/submits/$subfile"

#
# Submission
#

if $autosubmit; then

    condor_submit $PIPELINE_HOME/submits/$subfile
    emit_status $?

fi
