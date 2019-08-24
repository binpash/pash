#!/usr/bin/env dgsh
#
# SYNOPSIS Web log reporting
# DESCRIPTION
# Creates a report for a fixed-size web log file read from the standard input.
# Demonstrates the combined use of multipipe blocks, writeval and readval
# to store and retrieve values, and functions in the scatter block.
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

# Output the top X elements of the input by the number of their occurrences
# X is the first argument
toplist()
{
	uniq -c | sort -rn | head -$1
	echo
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

export -f toplist
export -f header


if [ -z "${DGSH_DRAW_EXIT}" ]
then
cat <<EOF
			WWW server statistics
			=====================

Summary
-------
EOF
fi

tee |
{{
	# Number of accesses
	echo -n 'Number of accesses: '
	dgsh-readval -l -s nAccess

	# Number of transferred bytes
	awk '{s += $NF} END {print s}' |
	tee |
	{{
		echo -n 'Number of Gbytes transferred: '
		awk '{print $1 / 1024 / 1024 / 1024}'

		dgsh-writeval -s nXBytes
	}}

	echo -n 'Number of hosts: '
	dgsh-readval -l -q -s nHosts

	echo -n 'Number of domains: '
	dgsh-readval -l -q -s nDomains

	echo -n 'Number of top level domains: '
	dgsh-readval -l -q -s nTLDs

	echo -n 'Number of different pages: '
	dgsh-readval -l -q -s nUniqPages

	echo -n 'Accesses per day: '
	dgsh-readval -l -q -s nDayAccess

	echo -n 'MBytes per day: '
	dgsh-readval -l -q -s nDayMB

	# Number of log file bytes
	echo -n 'MBytes log file size: '
	wc -c |
	awk '{print $1 / 1024 / 1024}'

	# Host names
	awk '{print $1}' |
	tee |
	{{
		# Number of accesses
		wc -l | dgsh-writeval -s nAccess

		# Sorted hosts
		sort |
		tee |
		{{

			# Unique hosts
			uniq |
			tee |
			{{
				# Number of hosts
				wc -l | dgsh-writeval -s nHosts

				# Number of TLDs
				awk -F. '$NF !~ /[0-9]/ {print $NF}' |
				sort -u |
				wc -l |
				dgsh-writeval -s nTLDs
			}}

			# Top 10 hosts
			{{
				 call 'header "Top 10 Hosts"'
				 call 'toplist 10'
			}}
		}}

		# Top 20 TLDs
		{{
			call 'header "Top 20 Level Domain Accesses"'
			awk -F. '$NF !~ /^[0-9]/ {print $NF}' |
			sort |
			call 'toplist 20'
		}}

		# Domains
		awk -F. 'BEGIN {OFS = "."}
		            $NF !~ /^[0-9]/ {$1 = ""; print}' |
		sort |
		tee |
		{{
			# Number of domains
			uniq |
			wc -l |
			dgsh-writeval -s nDomains

			# Top 10 domains
			{{
				 call 'header "Top 10 Domains"'
				 call 'toplist 10'
			}}
		}}
	}}

	# Hosts by volume
	{{
		call 'header "Top 10 Hosts by Transfer"'
		awk '    {bytes[$1] += $NF}
		END {for (h in bytes) print bytes[h], h}' |
		sort -rn |
		head -10
	}}

	# Sorted page name requests
	awk '{print $7}' |
	sort |
	tee |
	{{

		# Top 20 area requests (input is already sorted)
		{{
			 call 'header "Top 20 Area Requests"'
			 awk -F/ '{print $2}' |
			 call 'toplist 20'
		}}

		# Number of different pages
		uniq |
		wc -l |
		dgsh-writeval -s nUniqPages

		# Top 20 requests
		{{
			 call 'header "Top 20 Requests"'
			 call 'toplist 20'
		}}
	}}

	# Access time: dd/mmm/yyyy:hh:mm:ss
	awk '{print substr($4, 2)}' |
	tee |
	{{

		# Just dates
		awk -F: '{print $1}' |
		tee |
		{{

			# Number of days
			uniq |
			wc -l |
			tee |
			{{
				awk '
					BEGIN {
					"dgsh-readval -l -x -s nAccess" | getline NACCESS;}
					{print NACCESS / $1}' |
				dgsh-writeval -s nDayAccess

				awk '
					BEGIN {
					"dgsh-readval -l -x -q -s nXBytes" | getline NXBYTES;}
					{print NXBYTES / $1 / 1024 / 1024}' |
				dgsh-writeval -s nDayMB
			}}

			{{
				 call 'header "Accesses by Date"'
				 uniq -c
			}}

			# Accesses by day of week
			{{
				 call 'header "Accesses by Day of Week"'
				 sed 's|/|-|g' |
				 call '(date -f - +%a 2>/dev/null || gdate -f - +%a)' |
				 sort |
				 uniq -c |
				 sort -rn
			}}
		}}

		# Hour
		{{
			call 'header "Accesses by Local Hour"'
			awk -F: '{print $2}' |
			sort |
			uniq -c
		}}
	}}
	dgsh-readval -q -s nAccess
}} |
cat

