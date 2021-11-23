#!/bin/bash
# create bam files with regions
################### 1KG SAMPLES

IN=/home/dkarnikis/pash/evaluation/benchmarks/bio
IN_NAME=input/100G.txt #input.txt
OUT=/home/tammam1/pash/evaluation/benchmarks/bio/output

process_bio_s_line() {
        s_line=$1
        echo $s_line

        IN=/home/dkarnikis/pash/evaluation/benchmarks/bio
        IN_NAME=input.txt
        OUT=/home/tammam1/pash/evaluation/benchmarks/bio/output

        sample=$2 #(echo $s_line |cut -d " " -f 2);
        pop=$1 #(echo $s_line |cut -f 1 -d " ");
        link=$3 #(echo $s_line |cut -f 3 -d " ");
        ### correcting labeling of chromosomes so that all are 1,2,3.. instead of chr1,chr2 or chromosome1 etc
        echo 'Processing Sample '${IN}/input/$sample' ';
        # uniform the chromosomes in the file due to inconsistencies
        downloads/samtools-1.7/samtools view -H "${IN}/input/$sample".bam | sed -e 's/SN:\([0-9XY]\)/SN:chr\1/' -e 's/SN:MT/SN:chrM/' \
                | downloads/samtools-1.7/samtools reheader - "${IN}/input/$sample".bam > "${OUT}/$sample"_corrected.bam ;
        # create bai file 
        downloads/samtools-1.7/samtools index -b "${OUT}/$sample"_corrected.bam ;
        ### Isolating each relevant chromosome based on Gen_locs
        cut -f 2 ./Gene_locs.txt |sort |uniq |while read chr;
        do  
                echo 'Isolating Chromosome '$chr' from sample '${OUT}/$sample',  ';
                downloads/samtools-1.7/samtools view -b "${OUT}/$sample"_corrected.bam chr"$chr" > "${OUT}/$pop"_"$sample"_"$chr".bam ;
                echo 'Indexing Sample '$pop'_'${OUT}/$sample' ';
                downloads/samtools-1.7/samtools index -b "${OUT}/$pop"_"$sample"_"$chr".bam;
        done;
}

export -f process_bio_s_line
cat ${IN}/${IN_NAME} | xargs -n3 -I {} bash -c 'process_bio_s_line {}' 
