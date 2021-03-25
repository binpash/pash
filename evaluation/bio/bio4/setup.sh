PW=${PASH_TOP}/evaluation/bio/bio4/input
mkdir -p $PW
mkdir -p ${PASH_TOP}/evaluation/bio/bio4/output
cat ./input.txt |while read s_line;
	do
    sample=$(echo $s_line |cut -d " " -f 2);
    pop=$(echo $s_line |cut -f 1 -d " ");
    link=$(echo $s_line |cut -f 3 -d " ");
	wget -O "$PW/$sample".bam  "$link"; ##this part can be adjusted maybe
done;

