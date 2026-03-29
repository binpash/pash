#!/bin/bash
# tag: resize image

pure_func() {
    convert -resize 70% "-" "-"
}
export -f pure_func

for item in $(seq 1 "$ENTRIES"); do
    for image in $(cat "$PASH_TOP/evaluation/benchmarks/analytics/jpg_list.txt"); do
        cat "${IN}${image}" | pure_func > "${OUT}${image}.${item}"
    done
done

echo 'done';
