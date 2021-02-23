# Windows-like DIR command for the current directory. Nothing that couldn't be
# done with ls -l | awk. Demonstrates use of wrapped commands with no input (df,
# echo).
mkfifo a b c d e f a1 b1 c1 d1 e1 f1
cat a | awk '!/^total/ {print $6, $7, $8, $1, sprintf("%8d", $5), $9}' | cat > a1 &
cat b | wc -l | tr -d \\n | cat >b1 &
cat c | echo -n ' File(s) ' | cat > c1 &
# Tally number of bytes
cat d | awk '{s += $5} END {printf("%d bytes\n", s)}' | cat > d1 &
# Count number of directories
cat e | grep -c '^d' | tr -d \\n | cat > e1 &
# Print label for number of dirs and calculate free bytes
cat f | df -h . | awk '!/Use%/{print " Dir(s) " $4 " bytes free"}' | cat > f1 &
ls -n | tee a b c d e f> /dev/null &
sleep 1
cat a1 b1 c1 d1 e1 f1 
rm a1 b1 c1 d1 e1 f1 a b c d e f
