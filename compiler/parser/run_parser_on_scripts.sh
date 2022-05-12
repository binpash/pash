#! /bin/bash

SCRIPTS_DIR="../scripts/"

for script in "$SCRIPTS_DIR"*.sh
do
    echo "Parsing $script..."
    output=${script/"scripts"/"scripts/json"}.json
    ./parse_to_json.native "$script" > "$output"
done
