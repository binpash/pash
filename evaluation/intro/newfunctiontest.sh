cd "$(dirname $0)"

[ -z $PASH_TOP ] && {
  echo "PASH_TOP not set, maybe $(git rev-parse --show-toplevel)?"
  exit
}
FILE="input/100M.txt"

#test function to make a custom command for
our_func() {
        grep "$1"
}


cat "$FILE" | our_func choleric