###Download bam files
# create bam files with regions
################### 1KG SAMPLES
PW=$PASH_TOP/evaluation/scripts/input/bio4
cat ./Inds.txt |while read s_line;
	do
	sample=$(echo $s_line |cut -d " " -f 2);
	pop=$(echo $s_line |cut -f 1 -d " ");
	link=$(echo $s_line |cut -f 3 -d " ");
	### correcting labeling of chromosomes so that all are 1,2,3.. instead of chr1,chr2 or chromosome1 etc
	echo 'Processing Sample '$PW/$sample' ';
  # uniform the chromosomes in the file due to inconsistencies
	samtools view -H "$PW/$sample".bam | sed -e 's/SN:\([0-9XY]\)/SN:chr\1/' -e 's/SN:MT/SN:chrM/' \
    | samtools reheader - "$PW/$sample".bam > "$PW/$sample"_corrected.bam ;
	# create bai file 
	samtools index -b "$PW/$sample"_corrected.bam ;
	### Isolating each relevant chromosome based on Gen_locs
	cut -f 2 ./Gene_locs.txt |sort |uniq |while read chr;
		do  
			echo 'Isolating Chromosome '$chr' from sample '$PW/$sample',  ';
			samtools view -b "$PW/$sample"_corrected.bam chr"$chr" > "$PW/$pop"_"$sample"_"$chr".bam ;
			echo 'Indexing Sample '$pop'_'$PW/$sample' ';
			samtools index -b "$PW/$pop"_"$sample"_"$chr".bam;
			sleep 2
		done;
	rm "$PW/$sample"_corrected.bam;
	rm "$PW/$sample"_corrected.bam.bai;
	rm "$PW/$sample".bam
	done;
exit 1
