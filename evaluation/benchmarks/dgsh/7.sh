#!/bin/sh
# Automatically generated file
# Source file example/web-log-report.sh
#!/usr/bin/env sgsh -s /bin/bash
#
# SYNOPSIS Web log reporting
# DESCRIPTION
# Creates a report for a fixed-size web log file read from the standard input.
# Demonstrates the combined use of stores and named streams,
# the use of shell group commands and functions in the scatter block, and
# the use of cat(1) as a way to sequentially combine multiple streams.
# Used to measure throughput increase achieved through parallelism.
#
#  Copyright 2013 Diomidis Spinellis
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.
#
# FIXME need dataset
# Output the top X elements of the input by the number of their occurrences
# X is the first argument
toplist()
{
	uniq -c | sort -rn | head -$1
}

# Output the argument as a section header
header()
{
	echo
	echo "$1"
	echo "$1" | sed 's/./-/g'
}

# Consistent sorting
export LC_ALL=C

(

	export SGDIR=/tmp/sg-$$.0

	rm -rf $SGDIR

	# Cleanup on exit or interrupt
	cleanup()
	{
		SIGNAL=$1
		[ $SIGNAL = EXIT ] || echo sgsh interrupted. Cleaning up... 1>&2

		# Stop key-value stores
		
		# Kill processes we have launched in the background
		kill $SGPID 2>/dev/null

		# Remove temporary directory
		rm -rf "$SGDIR"

		# Propagate real signals and exit with non-0
		if [ $SIGNAL != EXIT ]
		then
			trap - $SIGNAL EXIT
			kill -s $SIGNAL $$
		fi

		# Exit with the original exit value
		exit

	}

	for sig in HUP INT QUIT TERM EXIT
	do
		trap "cleanup $sig" $sig
	done

	mkdir $SGDIR
	cat <&3 3<&-  >$SGDIR/npi-0.0.0
ln $SGDIR/npi-0.0.0 $SGDIR/npi-0.1.0
ln $SGDIR/npi-0.0.0 $SGDIR/npi-0.2.0
ln $SGDIR/npi-0.0.0 $SGDIR/npi-0.3.0
ln $SGDIR/npi-0.0.0 $SGDIR/npi-0.4.0
ln $SGDIR/npi-0.0.0 $SGDIR/npi-0.5.0
nXBytes=$(  {   awk '{s += $NF} END {print s}'
} <$SGDIR/npi-0.0.0 )
nLogBytes=$(  {   wc -c
} <$SGDIR/npi-0.1.0 )
 {   awk '{print $1}'
} <$SGDIR/npi-0.2.0 >$SGDIR/npi-1.0.0
ln $SGDIR/npi-1.0.0 $SGDIR/npi-1.1.0
ln $SGDIR/npi-1.0.0 $SGDIR/npi-1.2.0
ln $SGDIR/npi-1.0.0 $SGDIR/npi-1.3.0
nAccess=$(  {  wc -l
} <$SGDIR/npi-1.0.0 )
 {  sort
} <$SGDIR/npi-1.1.0 >$SGDIR/npi-2.0.0
ln $SGDIR/npi-2.0.0 $SGDIR/npi-2.1.0
 {  uniq
} <$SGDIR/npi-2.0.0 >$SGDIR/npi-3.0.0
ln $SGDIR/npi-3.0.0 $SGDIR/npi-3.1.0
nHosts=$(  {  wc -l
} <$SGDIR/npi-3.0.0 )
nTLD=$(  {  awk -F. '$NF !~ /[0-9]/ {print $NF}' |
					sort -u | wc -l
} <$SGDIR/npi-3.1.0 )
 {  {
				header 'Top 10 Hosts'
				toplist 10
			}
} <$SGDIR/npi-2.1.0 >$SGDIR/npfo-top10HostsByN.0
 {  {
			header 'Top 20 Level Domain Accesses'
			awk -F. '$NF !~ /^[0-9]/ {print $NF}' |
			sort |
			toplist 20
		}
} <$SGDIR/npi-1.2.0 >$SGDIR/npfo-top20TLD.0
 {  awk -F. 'BEGIN {OFS = "."}
		            $NF !~ /^[0-9]/ {$1 = ""; print}' | sort
} <$SGDIR/npi-1.3.0 >$SGDIR/npi-4.0.0
ln $SGDIR/npi-4.0.0 $SGDIR/npi-4.1.0
nDomain=$(  {  uniq | wc -l
} <$SGDIR/npi-4.0.0 )
 {  {
				header 'Top 10 Domains'
				toplist 10
			}
} <$SGDIR/npi-4.1.0 >$SGDIR/npfo-top10Domain.0
 {  {
		header 'Top 10 Hosts by Transfer'
		awk '    {bytes[$1] += $NF}
		END {for (h in bytes) print bytes[h], h}' |
		sort -rn |
		head -10
	}
} <$SGDIR/npi-0.3.0 >$SGDIR/npfo-top10HostsByVol.0
 {  awk '{print $7}' | sort
} <$SGDIR/npi-0.4.0 >$SGDIR/npi-5.0.0
ln $SGDIR/npi-5.0.0 $SGDIR/npi-5.1.0
ln $SGDIR/npi-5.0.0 $SGDIR/npi-5.2.0
 {  {
			header 'Top 20 Area Requests'
			awk -F/ '{print $2}' |
			toplist 20
		}
} <$SGDIR/npi-5.0.0 >$SGDIR/npfo-top20Area.0
nPages=$(  {  uniq | wc -l
} <$SGDIR/npi-5.1.0 )
 {  {
			header 'Top 20 Requests'
			toplist 20
		}
} <$SGDIR/npi-5.2.0 >$SGDIR/npfo-top20Request.0
 {  awk '{print substr($4, 2)}'
} <$SGDIR/npi-0.5.0 >$SGDIR/npi-6.0.0
ln $SGDIR/npi-6.0.0 $SGDIR/npi-6.1.0
 {  awk -F: '{print $1}'
} <$SGDIR/npi-6.0.0 >$SGDIR/npi-7.0.0
ln $SGDIR/npi-7.0.0 $SGDIR/npi-7.1.0
ln $SGDIR/npi-7.0.0 $SGDIR/npi-7.2.0
nDays=$(  {  uniq | wc -l
} <$SGDIR/npi-7.0.0 )
 {  {
				header 'Accesses by Date'
				uniq -c
			}
} <$SGDIR/npi-7.1.0 >$SGDIR/npfo-accessByDate.0
 {  {
				header 'Accesses by Day of Week'
				sed 's|/|-|g' |
				(date -f - +%a 2>/dev/null || gdate -f - +%a) |
				sort |
				uniq -c |
				sort -rn
			}
} <$SGDIR/npi-7.2.0 >$SGDIR/npfo-accessByDoW.0
 {  {
			header 'Accesses by Local Hour'
			awk -F: '{print $2}' |
			sort |
			uniq -c
		}
} <$SGDIR/npi-6.1.0 >$SGDIR/npfo-accessByHour.0
 { 
		cat <<EOF
			WWW server statistics
			=====================

Summary
-------
Number of accesses: $(echo ${nAccess})
Number of Gbytes transferred: $(expr $(echo ${nXBytes}) / 1024 / 1024 / 1024)
Number of hosts: $(echo ${nHosts})
Number of domains: $(echo ${nDomain})
Number of top level domains: $(echo ${nTLD})
Number of different pages: $(echo ${nPages})
Accesses per day: $(expr $(echo ${nAccess}) / $(echo ${nDays}))
MBytes per day: $(expr $(echo ${nXBytes}) / $(echo ${nDays}) / 1024 / 1024)
MBytes log file size: $(expr $(echo ${nLogBytes}) / 1024 / 1024)

EOF

}</dev/null  >$SGDIR/npfo-summary.0

# Gather the results

	# Gather all results together sequentially into a single report
	./sgsh-tee -f -I -i $SGDIR/npfo-summary.0 -i $SGDIR/npfo-top20Request.0 -i $SGDIR/npfo-top20Area.0 -i $SGDIR/npfo-top10HostsByN.0 -i $SGDIR/npfo-top10HostsByVol.0 -i $SGDIR/npfo-top10Domain.0 -i $SGDIR/npfo-top20TLD.0 -i $SGDIR/npfo-accessByDoW.0 -i $SGDIR/npfo-accessByHour.0 -i $SGDIR/npfo-accessByDate.0

)  3<&0 
