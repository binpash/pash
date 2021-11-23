# red color
RED='\033[0;31m'
# reset the color
NC='\033[0m'

IN=${BIO4:-$PASH_TOP/evaluation/benchmarks/bio}
IN_NAME=${IN_N:-input_all.txt}
if [[ $1 == "-c" ]]; then
    rm -rf *.bam
    rm -rf *.sam
    rm -rf ../output
    exit
fi

PW=${PASH_TOP}/evaluation/benchmarks/bio/input
mkdir -p $PW
mkdir -p ${PASH_TOP}/evaluation/benchmarks/bio/output
# install dependencies
pkgs='samtools'
if ! dpkg -s $pkgs >/dev/null 2>&1; then
    sudo apt-get install $pkgs
    var=$(dpkg -s samtools | grep Version)
    if [[ ! $var == "Version: 1.7-1" ]]; then                                
        printf "${RED}Invalid Samtools Version\n"                            
        printf "Samtools Version: 1.7 (using htslib 1.7-2) IS required{NC}\n"
    fi                                                                       
fi

cat ${IN}/${IN_NAME} |while read s_line;
	do
    
    sample=$(echo $s_line |cut -d " " -f 2);
    if [[ ! -f $sample ]]; then
        pop=$(echo $s_line |cut -f 1 -d " ");
        link=$(echo $s_line |cut -f 3 -d " ");
        wget -O "$PW/$sample".bam  "$link"; ##this part can be adjusted maybe
    fi
done;
