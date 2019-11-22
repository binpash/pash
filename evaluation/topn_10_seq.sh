# Top-N (1000) terms
cat $IN $IN $IN $IN $IN $IN $IN $IN $IN $IN | tr -cs A-Za-z '\n' | tr A-Z a-z | sort | uniq -c | sort -rn | sed ${N}q
