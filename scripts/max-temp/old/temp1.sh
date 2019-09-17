# Get maximum temperature--this program results to a Hadoop program on the order of hundreds of lines of code.
# (It's exploiting the fact that average temperatures everywhere are above 0:-)
paste -d '' <(yes "$DL" | head -n 11) <(seq 2005 2015)
seq 5 | { xargs -n 1 | grep 5}
for year in {2005..2015};
do
  # 2005
  # -rwxr-xr-x  1 nv  staff   257B Sep 19 15:51 incr.sh*
  # -rwxr-xr-x  1 nv  staff   155B Sep 19 13:25 sine.gz
  # -rwxr-xr-x  1 nv  staff   848B Sep 19 15:53 sq.sh*
  # -rw-r--r--  1 nv  staff   113K Sep 19 15:51 temp.txt
  curl "ftp://ftp.ncdc.noaa.gov/pub/data/noaa/$year/" | grep gz | awk '{print $NF}' | xargs -n 1 sh -c 'curl "ftp://ftp.ncdc.noaa.gov/pub/data/noaa/2015/$0" > $0"'
  cat *2015.gz | gunzip |  awk '{ print substr($0, 89, 4)}' | grep -v 9999 | sort -rn | head
done
