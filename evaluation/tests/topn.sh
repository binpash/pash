# Top-N (1000) terms
cat $IN | tr -cs A-Za-z '\n' | tr A-Z a-z | sort | uniq -c | sort -rn | sed ${N}q
