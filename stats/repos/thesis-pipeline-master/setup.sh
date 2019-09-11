# Dependencies
#
# Java          1.8.0_112
# GATK          3.6+
# Picard Tools  2.6.0
# Samtools      1.3.1
# BWA           7.15
# Seqtk         1.2-r94

#
# References / Locations
#

printf "\nSetting Pipeline Variables...\n"

# Set File System Root Directory
root_dir=/home/exacloud/lustre1

# Set User Directory
user_dir=/home/users/cordier

# Set Pipeline Directory
pipeline_dir=$root_dir/MRD_aml/ForBackup/pipeline

# Set Reference Genome Directory/File
reference_dir=$root_dir/MRD_aml/ForBackup/pipeline/reference/hg19

# Set Output Directory
output_dir=$user_dir/io # Set to home directory off lustre FS for FTP access

# Set Modules Directory
modules_dir=$pipeline_dir/procedures

# Set Validation Data Directory
validation_dir=$pipeline_dir/validation

# Set Logging Directory
logging_dir=$pipeline_dir/logs

# Set Auto Submit Directory
submit_dir=$pipeline_dir/submits

#
# Scaffolding
#

printf "\nScaffolding Directories...\n"

# Scaffold Pipeline Directory & Change Directory
cd $root_dir && mkdir -p $pipeline_dir && cd $pipeline_dir && mkdir -p $logging_dir/truth && mkdir -p $logging_dir && mkdir -p $submit_dir/logs

# Pipeline Scaffolding
mkdir -p set{1..6}/model/{rcorrector,bfc,musket,decgpu,lighter,bayeshammer,blessec,bloocoo,karect,kgem,quorum,seecer,shorah,nomodel,norealign}/param/{default,custom}/{logs,post-align,pre-align,markdup,recal,merged,modeled}
mkdir -p set{1..6}/model/{rcorrector,bfc,musket,decgpu,lighter,bayeshammer,blessec,bloocoo,karect,kgem,quorum,seecer,shorah,nomodel,norealign}/param/{default,custom}/modeled/{pairs,indexes}
mkdir -p set{1..6}/model/{rcorrector,bfc,musket,decgpu,lighter,bayeshammer,blessec,bloocoo,karect,kgem,quorum,seecer,shorah,nomodel,norealign}/param/{default,custom}/{post-align,pre-align}/fastq
mkdir -p set{1..6}/model/{rcorrector,bfc,musket,decgpu,lighter,bayeshammer,blessec,bloocoo,karect,kgem,quorum,seecer,shorah,nomodel,norealign}/param/{default,custom}/recal/{bqsr,nobqsr}/logs/{contest,mutect2,mutect}
mkdir -p set{1..6}/model/{rcorrector,bfc,musket,decgpu,lighter,bayeshammer,blessec,bloocoo,karect,kgem,quorum,seecer,shorah,nomodel,norealign}/param/{default,custom}/recal/bqsr/logs/bqsr
mkdir -p set{1..6}/{fastq,downloaded,tmp}
mkdir -p set{1..6}/downloaded/{intervals,metrics,original,split}
mkdir -p set{1..6}/downloaded/split/sorted
mkdir -p set{1..6}/fastq/{fastqc,split}
mkdir -p set{1..6}/fastq/split/{unpaired,unmerged}

# Clean Up Old pipeline.config
if [ -f $pipeline_dir/core/pipeline.config ]; then
    rm $pipeline_dir/core/pipeline.config
fi

# Create pipeline.state Write File
touch $pipeline_dir/core/pipeline.state

# Create pipeline.config File & Insert Shbang
touch $pipeline_dir/core/pipeline.config && echo -e "#! /usr/bin/bash\n" >> $pipeline_dir/core/pipeline.config

while true; do
    echo
    IFS= read -p "Setup Development Data Directory ($pipeline_dir/dev)? " yn
    case $yn in
        [Yy]* )  
            # Scaffold Development Directory
            mkdir -p dev/model/{rcorrector,bfc,musket,decgpu,lighter,bayeshammer,blessec,bloocoo,karect,kgem,quorum,seecer,shorah,nomodel,norealign}/param/{default,custom}/{logs,post-align,pre-align,markdup,recal,merged,modeled}
            mkdir -p dev/model/{rcorrector,bfc,musket,decgpu,lighter,bayeshammer,blessec,bloocoo,karect,kgem,quorum,seecer,shorah,nomodel,norealign}/param/{default,custom}/modeled/{pairs,indexes}
            mkdir -p dev/model/{rcorrector,bfc,musket,decgpu,lighter,bayeshammer,blessec,bloocoo,karect,kgem,quorum,seecer,shorah,nomodel,norealign}/param/{default,custom}/{post-align,pre-align}/fastq
            mkdir -p dev/model/{rcorrector,bfc,musket,decgpu,lighter,bayeshammer,blessec,bloocoo,karect,kgem,quorum,seecer,shorah,nomodel,norealign}/param/{default,custom}/recal/{bqsr,nobqsr}/logs/{contest,mutect2}
            mkdir -p dev/model/{rcorrector,bfc,musket,decgpu,lighter,bayeshammer,blessec,bloocoo,karect,kgem,quorum,seecer,shorah,nomodel,norealign}/param/{default,custom}/recal/bqsr/logs/bqsr
            mkdir -p dev/{fastq,downloaded,tmp}
            mkdir -p dev/downloaded/{intervals,metrics,original,split}
            mkdir -p dev/downloaded/split/sorted
            mkdir -p dev/fastq/{fastqc,split}
            mkdir -p dev/fastq/split/{unpaired,unmerged}
            break
            ;;
        [Nn]* ) 
            break
            ;;
        * ) echo "Requires yes [y] or no [n] input.";;
    esac
done

#
# Write Exports to pipeline.config
#

while true; do
    echo
    IFS= read -p "Write default pipeline exports to pipeline.config? " yn
    case $yn in
        [Yy]* ) (

                printf "Exporting Directory Locations...\n"

                # Header for Bash Profile
                echo -e '\n#\n# Pipeline Exports\n#' >> $pipeline_dir/core/pipeline.config

                # Write Pipeline Home Directory to pipeline.config
                echo -e "\nexport PIPELINE_HOME=$pipeline_dir" >> $pipeline_dir/core/pipeline.config
                # Set Pipeline Home Directory
                export PIPELINE_HOME=$pipeline_dir

                # Write Pipeline Modules Directory to pipeline.config
                echo -e "\nexport PIPELINE_MODS=$modules_dir" >> $pipeline_dir/core/pipeline.config
                # Set Pipeline Modules Directory
                export PIPELINE_MODS=$modules_dir

                # Write Pipeline Output Directory to pipeline.config
                echo -e "\nexport PIPELINE_OUT=$output_dir" >> $pipeline_dir/core/pipeline.config
                # Set Pipeline Output Directory
                export PIPELINE_OUT=$output_dir
                
                # Write Pipeline Reference Directory to pipeline.config
                echo -e "\nexport PIPELINE_REF=$reference_dir" >> $pipeline_dir/core/pipeline.config
                # Set Pipeline Reference Directory
                export PIPELINE_REF=$reference_dir
                
                # Write Pipeline Logging Directory to pipeline.config
                echo -e "\nexport PIPELINE_LOG=$logging_dir" >> $pipeline_dir/core/pipeline.config
                # Set Pipeline Logging Directory
                export PIPELINE_LOG=$logging_dir

            )
            break
            ;;
        [Nn]* ) 
            break
            ;;
        * ) echo "Requires yes [y] or no [n] input.";;
    esac
done

#
# Write Logging Management to pipeline.config
#

while true; do
    echo
    IFS= read -p "Enforce TMUX session logging via pipeline.config? " yn
    case $yn in
        [Yy]* )  
            # Inject Logging Automation Script to ~/.bash_profile
            echo -e '\n# Pipeline Logging Management\nif [[ $TERM = "screen" ]] && [[ $(ps -p $PPID -o comm=) = "tmux" ]] && [[ $PWD/ = $PIPELINE_HOME/* ]]; then\n\t# Prompt User / Read Input\n\tcd $PIPELINE_HOME/procedures\n\tIFS= read -p "Enter job logging name: " name\n\tmkdir $PIPELINE_HOME/logs 2> /dev/null\n\tlogname="tmux.$name.$(date "+%d%m%Y%H%M%S").log"\n\tscript -f $PIPELINE_LOG/${logname}\n\texit\nfi\n' >> $pipeline_dir/core/pipeline.config
            break
            ;;
        [Nn]* ) 
            break
            ;;
        * ) echo "Requires yes [y] or no [n] input.";;
    esac
done

#
# Inject Pipeline Config Import to ~/.bash_profile & ~/.bashrc
#

while true; do
    echo
    IFS= read -p "Set pipeline.config Reference in ~/.bash_profile & ~/.bashrc? " yn
    case $yn in
        [Yy]* )  
            # Inject pipeline.config to ~/.bash_profile & ~/.bashrc
            echo -e '\n#\n# Pipeline Config Reference\n#\n\n# Pipeline Config Link\nif [ -f $PIPELINE_HOME/core/pipeline.config ]; then\n\tsource $PIPELINE_HOME/core/pipeline.config\nfi' >> ~/.bash_profile
            echo -e '\n#\n# Pipeline Config Reference\n#\n\n# Pipeline Config Link\nif [ -f $PIPELINE_HOME/core/pipeline.config ]; then\n\tsource $PIPELINE_HOME/core/pipeline.config\nfi' >> ~/.bashrc
            break
            ;;
        [Nn]* ) 
            break
            ;;
        * ) echo "Requires yes [y] or no [n] input.";;
    esac
done

# Inject Import Script into pipeline.config
echo -e "#! /usr/bin/bash\n\n#\n# Pipeline Core Import\n#\n\nif [ -f $PIPELINE_HOME/core/pipeline.core ]; then\n\tsource $PIPELINE_HOME/core/pipeline.core\nfi\n" >> $pipeline_dir/core/pipeline.config

