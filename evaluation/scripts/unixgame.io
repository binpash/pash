# scripts from https://unixgame.io/
# https://github.com/psinghbh/softsec.github.io
cat input.txt | cut -d ' ' -f 2 | sort
cat input.txt | head -n 2 | cut -d ' ' -f 2
cat input.txt | cut -d ' ' -f 1 | sort | uniq -c | sort -r
cat input.txt | cut -d ' ' -f 4 | tr -d ','
cat input.txt | cut -d ' ' -f 2 | cut -c 1-1 | tr -d '\n' | tr '[A-Z]' '[a-z]'
cat input.txt | fmt -w1 | grep '\.' | wc -l
cat input.txt | fmt -w1 | grep 'x' | grep '\.' | wc -l
cat input.txt | fmt -w1 | grep 'x' | grep '\.' | cut -d '.' -f 2 | grep -v '[KQRBN]' | wc -l
cat input.txt | fmt -w1 | grep 'x' | grep '\.' | cut -d '.' -f 2 | grep '[KQRBN]' | cut -c 1-1 | sort | uniq -c | sort -nr
cat input.txt | fmt -w1 | grep 'x' | grep '\.' | cut -d '.' -f 2 | cut -c 1-1 | tr '[a-z]' 'P' | sort | uniq -c | sort -nr
cat input.txt | fmt -w1 | grep '\.' | cut -d '.' -f 2 | cut -c 1-1 | tr '[a-z]' 'P' | sort | uniq | tail -n 3 | sed 2d | sed 2d
cat input.txt | grep 'print' | cut -d '\"' -f 2 | cut -c 1-12
cat input.txt | awk '{print \$2, \$0}' | sort -nr | cut -d ' ' -f 2
cat input.txt | cut -f 1 | grep 'AT&T' | wc -l
cat input.txt | cut -f 2 | sort -n | uniq -c | sort -nr | head -n 1 | fmt -w1 | tail -n 1
cat input.txt | cut -f 4 | sort -n | cut -c 3-3 | uniq | sed s/\$/'0s'/
cat input.txt | fmt -w1 | grep 1969 | wc -l
cat input.txt | grep 'Bell' | awk 'length <= 45' | cut -d ',' -f 2
cat input.txt | grep '(' | cut -d '(' -f 2 | cut -d ')' -f 1 | head -n 1
cat input.txt | fmt -w1 | tr -c '[a-z][A-Z]' '\n' | sort | awk 'length >= 16'
cat input.txt | fmt -w1 | grep '[A-Z]' | tr '[a-z]' '\n' | grep '[A-Z]' | tr -d '\n' | cut -c 1-4
cat input.txt | cut -c 1-1 | tr -d '\n'
cat input.txt | cut -c 1-2 | tr -d '\n'
cat input.txt | fmt -w1 | grep '\"' | sed 4d | cut -d '\"' -f 2 | tr -d '\n'
cat input.txt | fmt -w1 | grep '[A-Z]' | sed 1d | sed 3d | sed 3d | tr '[a-z]' '\n' | grep '[A-Z]' | sed 3d | tr -c '[A-Z]' '\n' | tr -d '\n'
cat input.txt | sed 2d | sed 2d | fmt -w1 | grep '[A-Z]' | tr -c '[A-Z]' '\n' | tr -d '\n'
cat input.txt | tr -c '[a-z][A-Z]' '\n' | grep '[A-Z]' | sed 1d | sed 2d | sed 3d | sed 4d | tr -c '[A-Z]' '\n' | tr -d '\n'
cat input.txt | tr -c '[a-z][A-Z]' '\n' | grep '[A-Z]' | sed 1d | sed 1d | sed 2d | sed 3d | sed 5d | tr -c '[A-Z]' '\n' | tr -d '\n'
cat input.txt | sed 1d | grep 'Bell' | cut -f 2 | wc -l
cat input.txt | sed 1d | grep 'Bell' | cut -f 2
cat input.txt | grep 'Bell' | cut -f 2 | head -n 1 | fmt -w1 | cut -c 1-1 | tr -d '\n' | tr '[A-Z]' '[a-z]'
cat input.txt | grep 'UNIX' | cut -f 1
cat input.txt | cut -f 2 | cut -d ' ' -f 1 | sort | uniq -c | sort -nr | head -n 1 | fmt -w1 | sed 1d
