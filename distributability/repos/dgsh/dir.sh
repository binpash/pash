#!/usr/bin/env dgsh
#
# SYNOPSIS Directory listing
# DESCRIPTION
# Windows-like DIR command for the current directory.
# Nothing that couldn't be done with <code>ls -l | awk</code>.
# Demonstrates use of wrapped commands with no input (df, echo).
#
#  Copyright 2012-2017 Diomidis Spinellis
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

ls -n |
tee |
{{
	# Reorder fields in DIR-like way
	awk '!/^total/ {print $6, $7, $8, $1, sprintf("%8d", $5), $9}'

	# Count number of files
	wc -l | tr -d \\n

	# Print label for number of files
	echo -n ' File(s) '

	# Tally number of bytes
	awk '{s += $5} END {printf("%d bytes\n", s)}'

	# Count number of directories
	grep -c '^d' | tr -d \\n

	# Print label for number of dirs and calculate free bytes
	df -h . | awk '!/Use%/{print " Dir(s) " $4 " bytes free"}'
}} |
cat
