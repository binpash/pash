#!/bin/bash

# Consistent sorting
export LC_ALL=C

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

# Print initial header only if DGSH_DRAW_EXIT is not set
if [ -z "${DGSH_DRAW_EXIT}" ]
then
    cat <<EOF
			WWW server statistics
			=====================

Summary
-------
EOF
fi

## Initialize temporary files
file_initial=$(mktemp)
file_bytes=$(mktemp)
file_hosts=$(mktemp)
file_requests=$(mktemp)
file_sorted_hosts=$(mktemp)
file_unique_hosts=$(mktemp)
file_domains=$(mktemp)
file_day_count=$(mktemp)
file_dates=$(mktemp)
file_times=$(mktemp)

# This file will capture a large portion of the processed data to be reused in subsequent parts
cat > "$file_initial"

# Number of accesses
echo -n 'Number of accesses: '
wc -l "$file_initial"

# Total transferred bytes
awk '{s += $NF} END {print s}' "$file_initial" > "$file_bytes"
echo -n 'Number of Gbytes transferred: '
awk '{print $1 / 1024 / 1024 / 1024}' "$file_bytes"

# Process Host names
awk '{print $1}' "$file_initial" > "$file_hosts"

# Number of accesses
echo -n 'Number of accesses: '
wc -l < "$file_hosts"

# Sorted hosts
sort "$file_hosts" > "$file_sorted_hosts"

# Unique hosts
uniq "$file_sorted_hosts" > "$file_unique_hosts"
echo -n 'Number of hosts: '
wc -l < "$file_unique_hosts"

# Number of TLDs
echo -n 'Number of top level domains: '
awk -F. '$NF !~ /[0-9]/ {print $NF}' "$file_unique_hosts" | sort -u | wc -l


# Top 10 hosts
toplist 10 < "$file_sorted_hosts"

uniq -c "$file_sorted_hosts" | sort -rn | head -10
echo

# Top 20 TLDs
header "Top 20 Top Level Domains"

awk -F. '$NF !~ /^[0-9]/ {print $NF}' "$file_sorted_hosts" | sort | toplist 20
echo

# Domains
awk -F. 'BEGIN {OFS = "."} $NF !~ /^[0-9]/ {$1 = ""; print}' "$file_sorted_hosts" | sort > "$file_domains"

# Number of domains
echo -n 'Number of domains: '
uniq "$file_domains" | wc -l

# Top 10 domains
header "Top 10 domains"
toplist 10 < "$file_domains"

# Hosts by volume
header Top 10 Hosts by Transfer
awk '    {bytes[$1] += $NF}
END {for (h in bytes) print bytes[h], h}' "$file_initial" | sort -rn | head -10

# Sorted page name requests
awk '{print $7}' "$file_initial" | sort > "$file_requests"

# Top 20 area requests (input is already sorted)
header "Top 20 area requests"
awk -F/ '{print $2}' "$file_requests" | toplist 20
# Number of different pages
echo -n 'Number of different pages: '
cat "$file_requests" | uniq | wc -l

# Top 20 requests
header "Top 20 requests"
toplist 20 < "$file_requests"

# Access time: dd/mmm/yyyy:hh:mm:ss
awk '{print substr($4, 2)}' "$file_initial" > "$file_times"

# Just dates
awk -F: '{print $1}' "$file_times" > "$file_dates"

# Number of days
echo -n 'Accesses per day: '
uniq "$file_dates" | wc -l > "$file_day_count"
awk '
BEGIN {
    getline NACCESS < "'"$file_initial"'"
}
{print NACCESS / $1}' "$file_day_count"

echo -n 'MBytes per day: '
awk '
BEGIN {
    getline NXBYTES < "'"$file_bytes"'"
}
{print NXBYTES / $1 / 1024 / 1024}' "$file_day_count"

header "Accesses by Date"
uniq -c < "$file_dates"

# Accesses by day of week
header "Accesses by Day of Week"
sed 's|/|-|g' "$file_dates" | date -f - +%a 2>/dev/null | sort | uniq -c | sort -rn

# Accesses by Local Hour
header "Accesses by Local Hour"
awk -F: '{print $2}' "$file_times" | sort | uniq -c
