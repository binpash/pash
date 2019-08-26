# Get maximum temperature--this program results to a Hadoop program on the order of hundreds of lines of code.
# (It's exploiting the fact that average temperatures everywhere are above 0:-)
paste -d '' <(yes "$DL" | head -n 11) <(seq 2005 2015)
seq 5 | { xargs -n 1 | grep 5}
for year in {2005..2015};
do
  curl "ftp://ftp.ncdc.noaa.gov/pub/data/noaa/$year/"
  cat urls.dir | grep gz | awk '{print $NF}' | xargs -n 1 sh -c 'curl "ftp://ftp.ncdc.noaa.gov/pub/data/noaa/2015/$0" > $0"'
  cat *2015.gz | gunzip |  awk '{ print substr($0, 89, 4)}' | grep -v 9999 | sort -rn | head
done
