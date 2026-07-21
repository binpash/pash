#!/bin/bash

FILE="$PASH_TOP/../../README.md"

"cat" "$FILE" | "tr" "A-Z" "a-z"
'cat' "$FILE" | c"a"t | 'tr' 'a-z' 'A-Z'
