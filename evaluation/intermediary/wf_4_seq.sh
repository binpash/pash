# Word frequencies:
cat $IN $IN $IN $IN | tr -cs A-Za-z'\n' | tr A-Z a-z | sort | uniq -c | sort -rn | sed ${1}q
