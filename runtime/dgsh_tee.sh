
input=${1?"ERROR: dgsh-tee: No input file given"}
output=${2?"ERROR: dgsh-tee: No output file given"}
args=${@:3}

# Set a default DISH_TOP in this directory if it doesn't exist
PASH_TOP=${PASH_TOP:-$(git rev-parse --show-toplevel)}

$PASH_TOP/runtime/dgsh-tee -i "$input" -o "$output" $args

