#!/usr/bin/env dgsh
#
# SYNOPSIS C code metrics
# DESCRIPTION
# Process a directory containing C source code, and produce a summary
# of various metrics.
# Demonstrates nesting, commands without input.
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

{{
	# C and header code
	find "$@" \( -name \*.c -or -name \*.h \) -type f -print0 |
	tee |
	{{

		# Average file name length
		# Convert to newline separation for counting
		echo -n 'FNAMELEN: '
		tr \\0 \\n |
		# Remove path
		sed 's|^.*/||' |
		# Maintain average
		awk '{s += length($1); n++} END {
			if (n>0)
				print s / n;
			else
				print 0; }'

		xargs -0 /bin/cat |
		tee |
		{{
			# Remove strings and comments
			sed 's/#/@/g;s/\\[\\"'\'']/@/g;s/"[^"]*"/""/g;'"s/'[^']*'/''/g" |
			cpp -P |
			tee |
			{{
				# Structure definitions
				echo -n 'NSTRUCT: '
				egrep -c 'struct[   ]*{|struct[   ]*[a-zA-Z_][a-zA-Z0-9_]*[       ]*{'
				#}} (match preceding openings)

				# Type definitions
				echo -n 'NTYPEDEF: '
				grep -cw typedef

				# Use of void
				echo -n 'NVOID: '
				grep -cw void

				# Use of gets
	  			echo -n 'NGETS: '
	  			grep -cw gets

				# Average identifier length
				echo -n 'IDLEN: '
				tr -cs 'A-Za-z0-9_' '\n' |
				sort -u |
				awk '/^[A-Za-z]/ { len += length($1); n++ } END {
					if (n>0)
						print len / n;
					else
						print 0; }'
			}}

			# Lines and characters
			echo -n 'CHLINESCHAR: '
			wc -lc |
			awk '{OFS=":"; print $1, $2}'

			# Non-comment characters (rounded thousands)
			# -traditional avoids expansion of tabs
			# We round it to avoid failing due to minor
			# differences between preprocessors in regression
			# testing
			echo -n 'NCCHAR: '
			sed 's/#/@/g' |
			cpp -traditional -P |
			wc -c |
			awk '{OFMT = "%.0f"; print $1/1000}'

			# Number of comments
			echo -n 'NCOMMENT: '
			egrep -c '/\*|//'

			# Occurences of the word Copyright
			echo -n 'NCOPYRIGHT: '
			grep -ci copyright
		}}
	}}

	# C files
	find "$@" -name \*.c -type f -print0 |
	tee |
	{{
		# Convert to newline separation for counting
		tr \\0 \\n |
		tee |
		{{
			# Number of C files
			echo -n 'NCFILE: '
			wc -l

			# Number of directories containing C files
			echo -n 'NCDIR: '
			sed 's,/[^/]*$,,;s,^.*/,,' |
			sort -u |
			wc -l
		}}

		# C code
		xargs -0 /bin/cat |
		tee |
		{{
			# Lines and characters
			echo -n 'CLINESCHAR: '
			wc -lc |
			awk '{OFS=":"; print $1, $2}'

			# C code without comments and strings
			sed 's/#/@/g;s/\\[\\"'\'']/@/g;s/"[^"]*"/""/g;'"s/'[^']*'/''/g" |
			cpp -P |
			tee |
			{{
				# Number of functions
				echo -n 'NFUNCTION: '
				grep -c '^{'

				# Number of gotos
				echo -n 'NGOTO: '
				grep -cw goto

				# Occurrences of the register keyword
				echo -n 'NREGISTER: '
				grep -cw register

				# Number of macro definitions
				echo -n 'NMACRO: '
				grep -c '@[   ]*define[   ][   ]*[a-zA-Z_][a-zA-Z0-9_]*('
				# Number of include directives
				echo -n 'NINCLUDE: '
				grep -c '@[   ]*include'

				# Number of constants
				echo -n 'NCONST: '
				grep -ohw '[0-9][x0-9][0-9a-f]*' | wc -l

			}}
		}}
	}}

	# Header files
	echo -n 'NHFILE: '
	find "$@" -name \*.h -type f |
	wc -l

}} |
# Gather and print the results
cat
