dirs=(
    "oneliners"
    "unix50"
    "nlp"
    "analytics-mts"
    "dependency_untangling"
    "max-temp"
)

files=(
    "distr.res"
    "par.res"
    "seq.res"
    "hadoopstreaming.res"
    "distr_no_du.res"
    "par_no_du.res"
)

if [[ $# -eq 0 ]]; then
    echo "Usage: ./save_results.sh <dir_to_save_to>"
    exit 1
fi

save_to_dir=$1

echo saving to $save_to_dir

mkdir -p $save_to_dir

for dir in ${dirs[@]}
do
    mkdir -p $save_to_dir/$dir
    for file in ${files[@]}
    do
        file_path="$dir/$file"
        if [[ -f $file_path ]]; then
            cp $file_path "$save_to_dir/$file_path"
        fi
    done
done