# Top-N (1000) terms
N=1000
cat ./input.txt | tr -cs A-Za-z '\n' | tr A-Z a-z | sort | uniq -c | sort -rn | sed ${N}q

