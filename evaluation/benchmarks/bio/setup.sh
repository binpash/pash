if [[ $1 == "-c" ]]; then
    rm -rf input
    rm -rf output
    exit
fi
PW=${PASH_TOP}/evaluation/benchmarks/bio/input
mkdir -p $PW
mkdir -p ${PASH_TOP}/evaluation/benchmarks/bio/output
apt-get install samtools
cat ./input.txt |while read s_line;
	do
    
    sample=$(echo $s_line |cut -d " " -f 2);
    if [[ ! -f $sample ]]; then
        pop=$(echo $s_line |cut -f 1 -d " ");
        link=$(echo $s_line |cut -f 3 -d " ");
        wget -O "$PW/$sample".bam  "$link"; ##this part can be adjusted maybe
    fi
done;

