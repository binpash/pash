#!/bin/sh
IN=${IN:-10000.txt}
rm -f error.log
pkg_count=0
MIR_BIN=mir-sa
run_tests() {
    $MIR_BIN $1 -p >/dev/null 2>>error.log
}
export -f run_tests

rm -rf logs
mkdir -p logs

for item in node_modules/*;
do
    
    pkg_count=$((pkg_count + 1));
    #bash -c 'mir $0' $item
    run_tests $item > logs/$pkg_count.log
done

echo "done"
