#!/bin/bash

# Extracting honest pipelines from all of GitHub's aliases
# N.b.: Noisy data..

QUERY=$(cat <<'EOF'
.headers on
.mode csv
SELECT * FROM alias
EOF
)

# https://zenodo.org/record/3778825#.X9YamKpKjRZ
curl 'https://zenodo.org/record/3778825/files/results.db?download=1' > results.db # 3.2GB
sqlite3 results.db <(echo $QUERY) | csvcut -c 4 pipelines.csv | awk '{$1=$1};1' | sort | uniq tee >(
    # Schwartzian transform
    awk -F'|' '{print NF,$0}' file | sort -nr | cut -d' ' -f2- > likely-longest-pipelines.txt
  ) >(
    tr '|' '\n' | awk '{$1=$1};1' | awk '{print $1}' | tr -cs 'A-Za-z' '\n' | sort | uniq -c | sort -rn > freq-commands.txt
  )
