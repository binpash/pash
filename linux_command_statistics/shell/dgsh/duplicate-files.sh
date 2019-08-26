#!/usr/bin/env dgsh
#
# SYNOPSIS Find duplicate files
# DESCRIPTION
# List the names of duplicate files in the specified directory.
# Demonstrates the combination of streams with a relational join.
#
#  Copyright 2012-2013 Diomidis Spinellis
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

# Create list of files
find "$@" -type f |

# Produce lines of the form
# MD5(filename)= 811bfd4b5974f39e986ddc037e1899e7
xargs openssl md5 |

# Convert each line into a "filename md5sum" pair
sed 's/^MD5(//;s/)= / /' |

# Sort by MD5 sum
sort -k2 |

tee |
{{

	# Print an MD5 sum for each file that appears more than once
	awk '{print $2}' | uniq -d

	# Promote the stream to gather it
	cat
}} |
# Join the repeated MD5 sums with the corresponding file names
# Join expects two inputs, second will come from scatter
# XXX make streaming input identifiers transparent to users
join -2 2 |

# Output same files on a single line
awk '
BEGIN {ORS=""}
$1 != prev && prev {print "\n"}
END {if (prev) print "\n"}
{if (prev) print " "; prev = $1; print $2}'
