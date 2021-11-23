# encrypt all files in a directory 

IN=${IN:-./pcap/pcaps}
OUT=${OUT:-./out}
rm -rf logs_enc
mkdir -p logs_enc
run_tests() {
        echo 'key' | openssl enc -aes-256-cbc -e -md sha512 -in $item -out $OUT/$(basename $1).enc --pass stdin 
}
export -f run_tests
pkg_count=0
for item in ${IN}/*;
do
    pkg_count=$((pkg_count + 1));
    run_tests $item  > logs_enc/${pkg_count}.log
done

echo "done"
