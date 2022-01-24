#!/bin/bash
# create bam files with regions
################### 1KG SAMPLES
IN=${IN:-$PASH_TOP/evaluation/benchmarks/dependency_untangling/input}
SAMTOOLS_BIN=${IN}/deps/samtools-1.7/samtools
OUT=${OUT:-$PASH_TOP/evaluation/benchmarks/dependency_untangling/input/output/bio}
LOGS=${OUT}/logs
IN_NAME=${IN}/bio/100G.txt
GENE_LOCS=${IN}/bio/Gene_locs.txt
mkdir -p ${LOGS}
run_tests() {
    s_line=$(echo $1 | tr '@' ' ')
    pop=$(echo $s_line |cut -f 1 -d " ");
    sample=$(echo $s_line |cut -d " " -f 2);
    link=$(echo $s_line |cut -f 3 -d " ");
    ### correcting labeling of chromosomes so that all are 1,2,3.. instead of chr1,chr2 or chromosome1 etc
    echo 'Processing Sample '${IN}/bio/$sample' ';
    # uniform the chromosomes in the file due to inconsistencies
    $SAMTOOLS_BIN view -H "${IN}/bio/$sample".bam | sed -e 's/SN:\([0-9XY]\)/SN:chr\1/' -e 's/SN:MT/SN:chrM/' \
        | $SAMTOOLS_BIN reheader - "${IN}/bio/$sample".bam > "${OUT}/$sample"_corrected.bam  2> /dev/null
    # create bai file 
    $SAMTOOLS_BIN index -b "${OUT}/$sample"_corrected.bam 2> /dev/null
    ### Isolating each relevant chromosome based on Gen_locs
    cut -f 2 ${IN}/bio/Gene_locs.txt |sort |uniq |while read chr;
    do  
        echo 'Isolating Chromosome '$chr' from sample '${OUT}/$sample',  ';
        $SAMTOOLS_BIN view -b "${OUT}/$sample"_corrected.bam chr"$chr" > "${OUT}/$pop"_"$sample"_"$chr".bam 2> /dev/null
        echo 'Indexing Sample '$pop'_'${OUT}/$sample' ';
        $SAMTOOLS_BIN index -b "${OUT}/$pop"_"$sample"_"$chr".bam 2> /dev/null
    done;
}

export -f run_tests
data=$(cat ${IN_NAME} | tr ' ' '@')
pkg_count=0
for item in $data;
do
    pkg_count=$((pkg_count + 1));
    run_tests $item > "${LOGS}"/"${pkg_count}.log"
done

echo 'done';
