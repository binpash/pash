# create bam files with regions
################### 1KG SAMPLES
IN=${INPUT:-$PASH_TOP/evaluation/benchmarks/bio}
IN_NAME=${IN_N:-input_all.txt}
OUT=${OUTPUT:-$PASH_TOP/evaluation/benchmarks/bio/output}
cat ${IN}/${IN_NAME}|while read s_line;
  do
    sample=$(echo $s_line |cut -d " " -f 2);
    pop=$(echo $s_line |cut -f 1 -d " ");
    link=$(echo $s_line |cut -f 3 -d " ");
    ### correcting labeling of chromosomes so that all are 1,2,3.. instead of chr1,chr2 or chromosome1 etc
    echo 'Processing Sample '${IN}/input/$sample' ';
    # uniform the chromosomes in the file due to inconsistencies
    samtools view -H "${IN}/input/$sample".bam | sed -e 's/SN:\([0-9XY]\)/SN:chr\1/' -e 's/SN:MT/SN:chrM/' \
      | samtools reheader - "${IN}/input/$sample".bam > "${OUT}/$sample"_corrected.bam ;
    # create bai file 
    samtools index -b "${OUT}/$sample"_corrected.bam ;
    ### Isolating each relevant chromosome based on Gen_locs
    cut -f 2 ./Gene_locs.txt |sort |uniq |while read chr;
  do  
    echo 'Isolating Chromosome '$chr' from sample '${OUT}/$sample',  ';
    samtools view -b "${OUT}/$sample"_corrected.bam chr"$chr" > "${OUT}/$pop"_"$sample"_"$chr".bam ;
    echo 'Indexing Sample '$pop'_'${OUT}/$sample' ';
    samtools index -b "${OUT}/$pop"_"$sample"_"$chr".bam;
    #sleep 2
  done;
  #rm "${OUT}/$sample"_corrected.bam;
  #rm "${OUT}/$sample"_corrected.bam.bai;
  #rm "${OUT}/$sample".bam
done;
