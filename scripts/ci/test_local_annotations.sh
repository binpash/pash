#!/bin/bash

# Get repo root from PASH_TOP (which points to src/pash)
REPO_ROOT="${PASH_TOP}/../.."

# Set the annotation file path, this assumes annotations is a sibling repository
FILE="$REPO_ROOT/annotations/pash_annotations/annotation_generation/AnnotationGeneration.py"
LOG_FILE="$REPO_ROOT/pash.log"

# Ensure the annotation file exists
if [[ ! -f "$FILE" ]]; then
    echo "Error: Annotation file not found at $FILE"
    exit 1
fi

echo "Annotation file before test:"
grep '"cat": "Cat"' "$FILE" || echo "cat annotation not found (possibly already commented)"

# -------------------------------
# First run (annotation present)
# -------------------------------
echo "Running with cat annotation enabled..."
"$REPO_ROOT/pa.sh" --local-annotations-dir "$REPO_ROOT/annotations/" \
    --assert_all_regions_parallelizable \
    -c "cat /usr/share/dict/words | grep '^un' | wc -l"
first_run_status=$?
echo "First run exit code: $first_run_status"

if [[ "$first_run_status" -ne 0 ]]; then
    echo " Error: First run failed when annotation was present"
    exit 1
fi

# -----------------------------------
# Comment out "cat" annotation
# -----------------------------------
echo "Commenting out cat annotation in: $FILE"

##store original file for recovery later
cp "$FILE" "AnnotationGenerationCopy.py"

sed -i 's/^\(\s*\)"cat": "Cat",/\1# "cat": "Cat",/' "$FILE"

echo "Annotation file after commenting out:"
grep '"cat": "Cat"' "$FILE" || echo "cat annotation is now commented"

# ------------------------------------------
# Second run (annotation removed)
# ------------------------------------------
echo "Running with cat annotation removed..."
"$REPO_ROOT/pa.sh" --local-annotations-dir "$REPO_ROOT/annotations/" \
    --assert_all_regions_parallelizable \
    -c "cat /usr/share/dict/words | grep '^un' | wc -l"
second_run_status=$?
echo "Second run exit code: $second_run_status"

if [[ "$second_run_status" -eq 0 ]]; then
    echo "Error: Second run should have failed (no 'cat' annotation) but it succeeded"

    # Restore annotation before exiting
    sed -i 's/^\(\s*\)# "cat": "Cat",/\1"cat": "Cat",/' "$FILE"
    echo "Annotation restored after failure."

    exit 1
fi

# --------------------------------
# Restore original cat annotation
# --------------------------------
mv "AnnotationGenerationCopy.py" "$FILE"

echo "Annotation file after restoration:"
grep '"cat": "Cat"' "$FILE" || {
    echo "Error: Failed to restore annotation!"
    exit 1
}

echo "Test complete and cat annotation restored."

# Success
exit 0
