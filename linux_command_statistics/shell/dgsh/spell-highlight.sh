#!/usr/bin/env dgsh
#
# SYNOPSIS Highlight misspelled words
# DESCRIPTION
# Highlight the words that are misspelled in the command's first
# argument.
# Demonstrates stream processing with multipipes and
# the avoidance of pass-through constructs to avoid deadlocks.
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

export LC_ALL=C

tee |
{{
	# Find errors
	{{
		# Obtain list of words in text
		tr -cs A-Za-z \\n |
		tr A-Z a-z |
		sort -u

		# Ensure dictionary is compatibly sorted
		sort /usr/share/dict/words
	}} |
	# List errors as a set difference
	comm -23

	# Pass through text
	cat
}} |
grep --fixed-strings --file=- --ignore-case --color --word-regex --context=2
