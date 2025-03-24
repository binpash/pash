#!/bin/bash

# Define the file path
FILE="$PASH_TOP/../annotations/pash_annotations/annotation_generation/AnnotationGeneration.py"


# First run the script normally, make sure that the optimization occurred by checking that the FIFOs were made
$PASH_TOP/pa.sh -c "cat /usr/share/dict/words /usr/share/dict/words /usr/share/dict/words /usr/share/dict/words /usr/share/dict/words /usr/share/dict/words /usr/share/dict/words /usr/share/dict/words | grep '\(.\).*\1\(.\).*\2\(.\).*\3\(.\).*\4' | wc -l" -d 2 --local-annotations-dir ~/annotations/ --log_file "$PASH_TOP/pash.log"

#I initially tried to use --assert_compiler_success but that kept leading to the script henging when a certain section was not parrallelizable
if grep -Eq "rm_pash_fifos|mkfifo_pash_fifos" "$PASH_TOP/pash.log"; then
    echo "Success: the line was in the output when the cat annotations are available"
else
    echo "Error: The line is not found!"
    exit 1
fi  # ✅ Added missing `fi` here!

# Comment out the line containing '"cat": "Cat"' if not already commented
if [[ ! -f "$FILE" ]]; then
    echo "Error: File not found at $FILE"
    exit 1
fi

sed -i 's/^\(\s*\)"cat": "Cat",/\1# "cat": "Cat",/' "$FILE"
echo "Successfully commented out the 'cat' entry in $FILE"

# After removing the annotation for cat, the code should not be optimized, and no FIFOs will be made
$PASH_TOP/pa.sh -c "cat /usr/share/dict/words /usr/share/dict/words /usr/share/dict/words /usr/share/dict/words /usr/share/dict/words /usr/share/dict/words /usr/share/dict/words /usr/share/dict/words | grep '\(.\).*\1\(.\).*\2\(.\).*\3\(.\).*\4' | wc -l" -d 2 --local-annotations-dir ~/annotations/ --log_file "$PASH_TOP/pash.log"

if grep -Eq "rm_pash_fifos|mkfifo_pash_fifos" "$PASH_TOP/pash.log"; then
    echo "Error: The line was found!"
    exit 1
else
    echo "Success: The line is NOT in the output when the cat annotations are not available"
fi

# Restore the original "cat" annotation
sed -i 's/^\(\s*\)# "cat": "Cat",/\1"cat": "Cat",/' "$FILE"
