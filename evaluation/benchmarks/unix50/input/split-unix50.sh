#!/bin/bash

awk -v RS= '{print > (NR ".txt")}' unix50.sh

for file in *.txt; do
  fname=$(basename -- "$file")
  fscript="${fname%.*}".sh
  echo $fscript
  echo '#!/bin/bash' > $fscript

  echo 'export IN_PRE=${IN_PRE:-$PASH_TOP/evaluation/benchmarks/unix50/input}' >> $fscript
  input=$(grep -o 'IN..' $file)
  grep "^$(echo $input | xargs)=" unix50.sh >> $fscript
  cat $file >> $fscript
  echo '' >> $fscript
done

