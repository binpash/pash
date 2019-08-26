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

        # Access & Write Files With This Quality Model
        -q=*|--qualitymodel=*)
        qualitymodel="${i#*=}"
        shift
        ;;

        # Germline Tag
    	-N=*|--normal=*)
    	normal="${i#*=}"
    	shift
    	;;

        # Somatic Tag
    	-T=*|--tumor=*)
    	tumor="${i#*=}"
    	shift
    	;;

    # Optional Arguments With Defaults

        # n Reads Per GB or Memory
        -r=*|--reads=*)
        readsOpt="${i#*=}"
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

        # Estimate Cross Sample Contamination
        --contamination)
        contaminationOpt=true
        ;;

        # Trigger Debugging Available in Tools
        --debug)
        debugOpt=true
        ;;

        # Directory Cleanup (Voids All Other Parameters)
        --clean)
        cleanOpt=true
        ;;

    # Invalid Argument Handler

        *)
        printf "Invalid/Unused Parameter: $i"
        ;;
        
    esac
done

# Defaults if No Arguments Passed
readsDef="0"
ncoresDef="16"
memoryDef="6G"
contaminationDef=false
cleanDef=false
debugDef=false

# Set Optional Values
reads=${readsOpt:-$readsDef}
ncores=${ncoresOpt:-$ncoresDef}
memory=${memoryOpt:-$memoryDef}
contamination=${contaminationOpt:-$contaminationDef}
clean=${cleanOpt:-$cleanDef}
debug=${debugOpt:-$debugDef}

# Get Max Allowable Memory
allocMemory=${memory//[GgMmKk]/}
allocSize=${memory//[0-9]/}
maxMemory=$((allocMemory * ncores))$allocSize
 
# Max Reads in RAM - 200,000 per GB
maxReads=$((allocMemory * $reads))

# Use Default n Reads or User Defines
if [ "$maxReads" = "0" ]; then
    maxReads="default"
    readbuffersize=""
else
    readbuffersize="--read_buffer_size "$maxReads
fi

printf "\nPARAMETERS:
GATK Directory      = $GATK
Data File Prefix    = $fileprefix
Data Subset         = $subset
Experiment          = $experiment
Parameter Set       = $parameters
Recalibration Model = $qualitymodel
Contamination       = $contamination
Normal File Tag     = $normal
Tumor File Tag      = $tumor
Memory              = $memory
Cores               = $ncores
Max Memory          = $maxMemory
Max Reads in Memory = $maxReads
Debug               = $debug
Do Cleanup          = $clean
\n"

# Set Directories
recalDir=$PIPELINE_HOME/$subset/model/$experiment/param/$parameters/recal/$qualitymodel
tmpDir=$PIPELINE_HOME/$subset/tmp

# Tool Specific Debugging - GATK
monitorThreads=""
performanceLog=false
loggingLevel="INFO"

if $debug; then 
    monitorThreads="--monitorThreadEfficiency"
    performanceLog=true
    loggingLevel="DEBUG"
fi

format_status "Running Contamination Estimation & Mutect2 Script"

#
# If Contamination Not Specified, Run ContEst
#

if $contamination; then

    # State Check - Run Block if it Has Not Already Been Executed Successfully
    state="$fileprefix.$subset.$experiment.$parameters.$qualitymodel:MUTECT2:1"
    if !(has_state $state); then

        format_status "ContEst Start"
        # Define Command
        call="java -Xmx$maxMemory \
        -Djava.io.tmpdir=$tmpDir \
        -jar $GATK -T ContEst \
        --precision 0.001 \
        -R $PIPELINE_REF/Homo_sapiens_assembly19.fasta \
        -I:eval $recalDir/$fileprefix.$subset.$tumor.$experiment.$parameters.$qualitymodel.bam \
        -I:genotype $recalDir/$fileprefix.$subset.$normal.$experiment.$parameters.$qualitymodel.bam \
        -pf $PIPELINE_REF/hg19_population_stratified_af_hapmap_3.3.cleaned.vcf \
        -isr INTERSECTION \
        --population ALL \
        --log_to_file $recalDir/logs/contest/log_$experiment-cont_est_recal.txt \
        -o $recalDir/logs/contest/cont_est_recal_$experiment.txt \
        --logging_level $loggingLevel \
        $monitorThreads $readbuffersize" # Additional Optional Args
        # Print & Call
        format_status "Command:\n$call"
        eval $call


        # Update State on Exit
        put_state $? $state
        format_status "ContEst Complete"

    fi

    #
    # Mutect2
    #

    # State Check - Run Block if it Has Not Already Been Executed Successfully
    state="$fileprefix.$subset.$experiment.$parameters.$qualitymodel:MUTECT2:2"
    if !(has_state $state); then

        # Get Contamination Fraction
        contaminationPercent=$(awk -F '\t' 'NR >=2 {print $4}'  $recalDir/logs/contest/cont_est_recal_$experiment.txt)
        contamination=$(python -c "print($contaminationPercent/100.0)")
        format_status "Proportion Contamination: $contamination"

        format_status "MuTect2 Start"
        # Define Command
        call="java -Xmx$memory \
        -Djava.io.tmpdir=$tmpDir \
        -jar $GATK -T MuTect2 \
        -R $PIPELINE_REF/Homo_sapiens_assembly19.fasta \
        -I:tumor $recalDir/$fileprefix.$subset.$tumor.$experiment.$parameters.$qualitymodel.bam \
        -I:normal $recalDir/$fileprefix.$subset.$normal.$experiment.$parameters.$qualitymodel.bam \
        --dbsnp $PIPELINE_REF/dbsnp_138.hg19_modified.vcf \
        --cosmic $PIPELINE_REF/b37_cosmic_v54_120711_modified.vcf \
        --tumor_lod 10.0 \
        --contamination_fraction_to_filter $contamination \
        -o $recalDir/logs/mutect2/$fileprefix.$subset.$experiment.$parameters.$qualitymodel.raw.snps.indels.vcf \
        --log_to_file $recalDir/logs/mutect2/log_mutect2_$experiment.txt \
        --graphOutput $recalDir/logs/mutect2/assembly_graph_info.txt \
        -nct $ncores \
        --logging_level $loggingLevel \
        $monitorThreads $readbuffersize" # Additional Optional Args
        # Print & Call
        format_status "Command:\n$call"
        eval $call

        # Update State on Exit
        put_state $? $state
        format_status "MuTect2 Complete"

    fi

else

    #
    # Mutect2
    #

    # State Check - Run Block if it Has Not Already Been Executed Successfully
    state="$fileprefix.$subset.$experiment.$parameters.$qualitymodel:MUTECT2:2"
    if !(has_state $state); then

        format_status "MuTect2 Start"
        # Define Command
        call="java -Xmx$memory \
        -Djava.io.tmpdir=$tmpDir \
        -jar $GATK -T MuTect2 \
        -R $PIPELINE_REF/Homo_sapiens_assembly19.fasta \
        -I:tumor $recalDir/$fileprefix.$subset.$tumor.$experiment.$parameters.$qualitymodel.bam \
        -I:normal $recalDir/$fileprefix.$subset.$normal.$experiment.$parameters.$qualitymodel.bam \
        --dbsnp $PIPELINE_REF/dbsnp_138.hg19_modified.vcf \
        --cosmic $PIPELINE_REF/b37_cosmic_v54_120711_modified.vcf \
        --tumor_lod 10.0 \
        --contamination_fraction_to_filter 0.0 \
        -o $recalDir/logs/mutect2/$fileprefix.$subset.$experiment.$parameters.$qualitymodel.raw.snps.indels.vcf \
        --log_to_file $recalDir/logs/mutect2/log_mutect2_$experiment.txt \
        --graphOutput $recalDir/logs/mutect2/assembly_graph_info.txt \
        -nct $ncores \
        --logging_level $loggingLevel \
        $monitorThreads $readbuffersize" # Additional Optional Args
        # Print & Call
        format_status "Command:\n$call"
        eval $call

        # Update State on Exit
        put_state $? $state
        format_status "MuTect2 Complete"

    fi

fi

# 
# Copy VCFS to User I/O Directory
# 

# State Check - Run Block if it Has Not Already Been Executed Successfully
state="$fileprefix.$subset.$experiment.$parameters.$qualitymodel:MUTECT2:3"
if !(has_state $state); then

    format_status "Copying VCFs to I/O Directory"
    # Define Command
    call="cp $recalDir/logs/mutect2/$fileprefix.$subset.$experiment.$parameters.$qualitymodel.raw.snps.indels.vcf /home/users/$USER/io/"
    # Print & Call
    format_status "Command:\n$call"
    eval $call

    # Update State on Exit
    put_state $? $state
    format_status "VCFs Copied to I/O Directory"

fi

format_status "Done"

