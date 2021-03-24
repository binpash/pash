# count how many times each file type exist in a directory
# 1GB files take about 14s on i7 8th gen
get_type_count()
(
    find $1 -type f | while read f; do echo ""${f##*.}""; done | sed ""/^\s*$/d"" | sort | uniq -c | sort -rn  > $2
)

get_type_count input/tmp $PWD/output/get_type_count_res
