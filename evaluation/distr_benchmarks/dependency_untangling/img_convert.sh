#!/bin/bash
# tag: resize image 
IN=${JPG:-/dependency_untangling/jpg}
OUT=${OUT:-$PASH_TOP/evaluation/distr_benchmarks/dependency_untangling/input/output/jpg}
mkdir -p ${OUT}

pure_func () {
     convert -resize 70% "-" "-"
}
export -f pure_func

for i in $(hdfs dfs -ls -C ${IN}/*.jpg); 
do 
    out=$OUT/$(basename -- $i)
    hdfs dfs -cat $i | pure_func > $out; 
done

echo 'done';
