#!/bin/bash
# create bam files with regions
################### 1KG SAMPLES
set -e
IN=${IN:-$PASH_TOP/evaluation/benchmarks/for-loops/input}
SAMTOOLS_BIN=${IN}/samtools-1.7/samtools
OUT=${OUT:-$PASH_TOP/evaluation/benchmarks/for-loops/input/output/bio}
LOGS=${OUT}/logs
IN_NAME=${IN}/100G.txt
GENE_LOCS=${IN}/Gene_locs.txt
mkdir -p ${LOGS}
run_tests() {
    #s_line=$1
    pop=$1 #(echo $s_line |cut -f 1 -d " ");
    sample=$2 #(echo $s_line |cut -d " " -f 2);
    link=$3 #(echo $s_line |cut -f 3 -d " ");
    ### correcting labeling of chromosomes so that all are 1,2,3.. instead of chr1,chr2 or chromosome1 etc
    echo 'Processing Sample '${IN}/bio/$sample' ';
    # uniform the chromosomes in the file due to inconsistencies
    $SAMTOOLS_BIN view -H "${IN}/bio/$sample".bam | sed -e 's/SN:\([0-9XY]\)/SN:chr\1/' -e 's/SN:MT/SN:chrM/' \
        | $SAMTOOLS_BIN reheader - "${IN}/bio/$sample".bam > "${OUT}/$sample"_corrected.bam ;
    # create bai file 
    $SAMTOOLS_BIN index -b "${OUT}/$sample"_corrected.bam ;
    ### Isolating each relevant chromosome based on Gen_locs
    cut -f 2 ${IN}/Gene_locs.txt |sort |uniq |while read chr;
    do  
        echo 'Isolating Chromosome '$chr' from sample '${OUT}/$sample',  ';
        $SAMTOOLS_BIN view -b "${OUT}/$sample"_corrected.bam chr"$chr" > "${OUT}/$pop"_"$sample"_"$chr".bam ;
        echo 'Indexing Sample '$pop'_'${OUT}/$sample' ';
        $SAMTOOLS_BIN index -b "${OUT}/$pop"_"$sample"_"$chr".bam;
    done;
}

export -f run_tests
data=$(head -n2 ${IN_NAME} | tr ' ' '@')
for pkg in $data;
do
    arg0=$(echo $pkg | cut -d "@" -f 1);
    arg1=$(echo $pkg | cut -d "@" -f 2);
    arg2=$(echo $pkg | cut -d "@" -f 3);
    run_tests $arg0 $arg1 $arg2 > ${LOGS}/$arg0.log
done

echo 'done';
