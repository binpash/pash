# sdsh
Scalable Distributed Shell or something!

The key insight here is that the key shell abstractions (streams and pipes) are trivially distributed in a scalable and fault-tolerant manner---for example, similar to Spark's RDD's.

# Possible Homework Breakdown

* Create a few categories---trivially distributable (e.g., pure functions), might require some coordination, 

* Break down most primitives 

  * maybe some shells (such as Plan 9's `rc` shell) have primitives that are naturally more amenable to distribution than others.

* What about variable sharing (mostly write once / read many times)

* Other constructs such as `if` and `for`?

* Need to Redirect input / output streams

* We need to look into an extensibility. What if the user, in their environment, have access to primitives not "known" to the shell?

# A Few Examples

Simple pipeline examples:

Word Count
```bash
cat $INPUT | wc -w 
```

Grep -- search a substring within a string
```bash
cat $INPUT | grep '[a-zA-Z0-9]\+@[a-zA-Z0-9]\+\.[a-z]\{2,\}' 
```

Word frequencies:
```bash
cat $INPUT | tr -cs A-Za-z'\n' | tr A-Z a-z | sort | uniq -c | sort -rn | sed ${1}q
```

Top-N (1000) terms
```bash
cat $INPUT | tr -cs A-Za-z '\n' | tr A-Z a-z | sort | uniq -c | sort -rn | sed 1000q
```

N-Grams (for N==2, bigrams)
```bash
cat $INPUT | tr -cs A-Za-z '\n' | tr A-Z a-z > tokens.txt && tail +2 tokens.txt > next.txt && paste tokens.txt next.txt > bigrams.txt && cat bigrams.txt | sort | uniq > results
```

Get maximum temperature--this program results to a Hadoop program on the order of hundreds of lines of code.
(It's exploiting the fact that average temperatures everywhere are above 0:-)
```bash
paste -d '' <(yes "$DL" | head -n 11) <(seq 2005 2015)
seq 5 | { xargs -n 1 | grep 5}
for year in {2005..2015}
  curl "ftp://ftp.ncdc.noaa.gov/pub/data/noaa/$year/"
  cat urls.dir | grep gz | awk '{print $NF}' | xargs -n 1 sh -c 'curl "ftp://ftp.ncdc.noaa.gov/pub/data/noaa/2015/$0" > $0"'
  cat *2015.gz | gunzip |  awk '{ print substr($0, 89, 4)}' | grep -v 9999 | sort -rn | head
done
```

The same, written in a different way
```bash
curl $url/$year | grep gz | awk '{print $NF}' | { xargs -n 1 curl > $year }
for yr in {2005..2015}; do curl $url/$year | grep gz | awk '{print $NF}' | { xargs -n 1 | gunzip | awk '{print substr($0, 89, 4)}' | grep -v 9999 | sort -rn | head -n 1 | sed "s/^/$yr: /" } done
```

Fine markdown files, compile them, and expose them over the network
(showcases piece-wise comments)
```bash
find . '*.md' | # Parallelizable, given a distributed FS
    xargs mdc | # xargs is higher-order command; trivially parallelizable
    nc -l 80    # netcat could default-but-configurably parallelizable
```

# (Some) Related Work

In shell scripting:
* Attempts on distributed shells: `rc`, 
* Parallel commands: GNU's `parallel`, 
* Tons of parallel scripting: Lsf, Swift, ..

Outside shell scripting:
* Cluster computing: Hadoop, esp. Spark
* Stream processing:
