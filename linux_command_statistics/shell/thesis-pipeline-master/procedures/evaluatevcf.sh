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

        -f=*|--fileprefix=*)
        fileprefix="${i#*=}"
        shift # Access & Write Files With This Prefix
        ;;
        -s=*|--subset=*)
        subset="${i#*=}"
        shift # Access & Write Files With This Subset
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

        -m=*|--memory=*)
        memoryOpt="${i#*=}"
        shift # Per Core Memory Requirement
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
memoryDef="6G"

# Set Optional Values
memory=${memoryOpt:-$memoryDef}

# Get Max Allowable Memory
allocMemory=${memory//[GgMmKk]/}
allocSize=${memory//[0-9]/}
maxMemory=$memory

printf "\nPARAMETERS:
GATK Directory      = $GATK
Data File Prefix    = $fileprefix
Data Subset         = $subset
Experiment          = $experiment
Parameter Set       = $parameters
Recalibration Model = $qualitymodel
Memory              = $memory
\n"

# Set Directories
recalDir=$PIPELINE_HOME/$subset/model/$experiment/param/$parameters/recal/$qualitymodel
tmpDir=$PIPELINE_HOME/$subset/tmp

# Assign Alternate Quality Model
if [[ "$qualitymodel" == "no"* ]]; then
    altqualitymodel=${qualitymodel//[NonoNOnO]/}
else
    altqualitymodel="no$qualitymodel"
fi

#
# Variant Evaluation
#

format_status "Running Variant Evaluation"

# State Check - Run Block if it Has Not Already Been Executed Successfully
state="$fileprefix.$subset.$experiment.$parameters.$qualitymodel:EVALVCF:1"
if !(has_state $state); then

    format_status "ContEst Start"
    # Define Command
    call="java -Xmx$maxMemory \
    -Djava.io.tmpdir=$tmpDir \
    -jar $GATK -T VariantEval \
    -R $PIPELINE_REF/Homo_sapiens_assembly19.fasta \
    --eval:set1 $recalDir/logs/mutect2/$fileprefix.$subset.$experiment.$parameters.$qualitymodel.raw.snps.indels.vcf \
    --comp $recalDir/logs/mutect2/$fileprefix.$subset.$experiment.$parameters.$altqualitymodel.raw.snps.indels.vcf \
    --goldStandard $PIPELINE_HOME/validation/truth/$subset.truth.vcf" 
    # Print & Call
    format_status "Command:\n$call"
    eval $call

    # Update State on Exit
    put_state $? $state
    format_status "EvalVCF Complete"

fi

format_status "Done"

