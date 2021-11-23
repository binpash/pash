# compress all files in a directory

IN=${IN:-./pcap/pcaps}
OUT=${OUT:-./out}
mkdir -p out
run_tests() {
        tar -czf $OUT/$(basename $1).tar.gz $1
}

export -f run_tests
pkg_count=0
for item in ${IN}/*;
do
    pkg_count=$((pkg_count + 1));
    run_tests $item > logs/$pkg_count.log
done

echo "done"
