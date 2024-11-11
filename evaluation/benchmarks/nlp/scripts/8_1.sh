#!/bin/bash
# tag: sort_words_by_num_of_syllables
# set -e

mkdir -p "$OUT"

inputs="
61-0.txt
61.txt
655.txt
660.txt
665.txt
66600-0.txt
66601-0.txt
66602-0.txt
66603-0.txt
66604-0.txt
66605-0.txt
66606-0.txt
66607-0.txt
66608-0.txt
66609-0.txt
6660.txt
6661.txt
66620-0.txt
66621-0.txt
66622-0.txt
66623-0.txt
66624-0.txt
66625-0.txt
66626-0.txt
66627-0.txt
66628-0.txt
66629-0.txt
6662.txt
66630-0.txt
66631-0.txt
66632-0.txt
66633-0.txt
66634-0.txt
66635-0.txt
66636-0.txt
66637-0.txt
66638-0.txt
66639-0.txt
6663.txt
66640-0.txt
6664-0.txt
66641-0.txt
66642-0.txt
66643-0.txt
66645-0.txt
66646-0.txt
66647-0.txt
66648-0.txt
66649-0.txt
6664.txt
66650-0.txt
66651-0.txt
66652-0.txt
66653-0.txt
66654-0.txt
66655-0.txt
66656-0.txt
66657-0.txt
66658-0.txt
6665-8.txt
66659-0.txt
66660-0.txt
6666-0.txt
66661-0.txt
66662-0.txt
66663-0.txt
66664-0.txt
66665-0.txt
66666-0.txt
66667-0.txt
66668-0.txt
66669-0.txt
66670-0.txt
6667-0.txt
66671-0.txt
66672-0.txt
66673-0.txt
66674-0.txt
66675-0.txt
66676-0.txt
66677-0.txt
66678-0.txt
66679-0.txt
66680-0.txt
66681-0.txt
66682-0.txt
66683-0.txt
66684-0.txt
66685-0.txt
66686-0.txt
66687-0.txt
66688-0.txt
66689-0.txt
66690-0.txt
66691-0.txt
66692-0.txt
66693-0.txt
66694-0.txt
66695-0.txt
66696-0.txt
66697-0.txt
66698-0.txt
66699-0.txt
6669.txt
666.txt
fchld10.txt
gdlns10.txt
helb10.txt
lttlc10.txt
manif11.txt
manif12.txt
pgw050ab.txt
pgw050mo.txt
pgw050pq.txt
pgwab04.txt
pgwmo04.txt
pgwpq04.txt
ppow10.txt
teop210.txt
wvr1210.txt
"

pure_func() {
    input=$1
    TEMPDIR=$(mktemp -d)
    cat > ${TEMPDIR}/${input}.words
    tr -sc '[AEIOUaeiou\012]' ' ' < ${TEMPDIR}/${input}.words | awk '{print NF}' > ${TEMPDIR}/${input}.syl
    paste ${TEMPDIR}/${input}.syl ${TEMPDIR}/${input}.words | sort -nr | sed 5q
    rm -rf ${TEMPDIR}
}
export -f pure_func
for input in $(echo $inputs | tr " " "\n" | head -n ${ENTRIES})
do
    python3 "$SERVERLESS_RUNTIME_DIR/s3-get-object.py" "$IN/$input" /dev/stdout |
    cat  | tr -c 'A-Za-z' '[\n*]' | grep -v "^\s*$" | sort -u | pure_func $input > /dev/stdout |
    python3 "$SERVERLESS_RUNTIME_DIR/s3-put-object.py" "$OUT/$input.out" /dev/stdin dummy # dummy object_id, not used for bash
done

echo 'done';
# rm -rf "$OUT"
