# works on pash
# Adapted from the DGSH
# https://github.com/dspinellis/dgsh/blob/master/example/reorder-columns.sh

mkfifo a b c d
cat a  | cut -d , -f 5-6 - | cat > c & 
cat b  | cut -d , -f 2-4 - | cat > d &
cat $1 | tee a b > /dev/null  &
paste -d , c d
rm a b c d
