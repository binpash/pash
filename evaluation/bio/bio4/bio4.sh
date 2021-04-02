###Download bam files
# create bam files with regions
################### 1KG SAMPLES
INPUT=${INPUT:-$PASH_TOP/evaluation/bio/bio4/input}
OUTPUT=${OUTPUT:-$PASH_TOP/evaluation/bio/bio4/output}

cat ${INPUT}/../input.txt |while read s_line;
	do
	sample=$(echo $s_line |cut -d " " -f 2);
	pop=$(echo $s_line |cut -f 1 -d " ");
	link=$(echo $s_line |cut -f 3 -d " ");
	### correcting labeling of chromosomes so that all are 1,2,3.. instead of chr1,chr2 or chromosome1 etc
	echo 'Processing Sample '${INPUT}/$sample' ';
  # uniform the chromosomes in the file due to inconsistencies
	samtools view -H "${INPUT}/$sample".bam | sed -e 's/SN:\([0-9XY]\)/SN:chr\1/' -e 's/SN:MT/SN:chrM/' \
    | samtools reheader - "${INPUT}/$sample".bam > "${OUTPUT}/$sample"_corrected.bam ;
	# create bai file 
	samtools index -b "${OUTPUT}/$sample"_corrected.bam ;
	### Isolating each relevant chromosome based on Gen_locs
	cut -f 2 ./Gene_locs.txt |sort |uniq |while read chr;
		do  
			echo 'Isolating Chromosome '$chr' from sample '${OUTPUT}/$sample',  ';
			samtools view -b "${OUTPUT}/$sample"_corrected.bam chr"$chr" > "${OUTPUT}/$pop"_"$sample"_"$chr".bam ;
			echo 'Indexing Sample '$pop'_'${OUTPUT}/$sample' ';
			samtools index -b "${OUTPUT}/$pop"_"$sample"_"$chr".bam;
			sleep 2
		done;
	rm "${OUTPUT}/$sample"_corrected.bam;
	rm "${OUTPUT}/$sample"_corrected.bam.bai;
	#rm "${OUTPUT}/$sample".bam
	done;
exit 1
