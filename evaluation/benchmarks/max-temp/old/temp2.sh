# The same as `temp1.sh`, written in a different way
curl $url/$year | grep gz | awk '{print $NF}' | { xargs -n 1 curl > $year }
for yr in {2005..2015}; do curl $url/$year | grep gz | awk '{print $NF}' | { xargs -n 1 | gunzip | awk '{print substr($0, 89, 4)}' | grep -v 9999 | sort -rn | head -n 1 | sed "s/^/$yr: /" } done
