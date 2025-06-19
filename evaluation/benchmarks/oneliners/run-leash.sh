
set -x

# time IN=oneliners/inputs/1M.txt OUT=oneliners/outputs/nfa-regex-1M-leash-w256-lambda-only- $PASH_TOP/pa.sh --serverless_exec -w256 scripts/nfa-regex.sh
# time IN=oneliners/inputs/1G.txt OUT=oneliners/outputs/nfa-regex-1G-leash-w256-lambda-only- $PASH_TOP/pa.sh --serverless_exec -w256 scripts/nfa-regex.sh
time IN=oneliners/inputs/3G.txt OUT=oneliners/outputs/nfa-regex-3G-leash-w256-lambda-only- $PASH_TOP/pa.sh --serverless_exec -w256 scripts/nfa-regex.sh

