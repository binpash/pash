#!/bin/bash

## Initialize the necessary temporary files
file1=$(mktemp)
file2=$(mktemp)
file3=$(mktemp)

# Create list of files
find "$@" -type f |

# Produce lines of the form
# MD5(filename)= 811bfd4b5974f39e986ddc037e1899e7
xargs openssl md5 |

# Convert each line into a "filename md5sum" pair
sed 's/^MD5(//;s/)= / /' |

# Sort by MD5 sum
sort -k2 > "$file1"

# Print an MD5 sum for each file that appears more than once
awk '{print $2}' < "$file1" | uniq -d > "$file2"


# Join the repeated MD5 sums with the corresponding file names
# Join expects two inputs, second will come from scatter
# XXX make streaming input identifiers transparent to users
join -2 2 "$file2" "$file1" |

# Output same files on a single line
awk '
BEGIN {ORS=""}
$1 != prev && prev {print "\n"}
END {if (prev) print "\n"}
{if (prev) print " "; prev = $1; print $2}'
