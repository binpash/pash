paste -d '+' <(cat * | wc | tr -s ' '  '\n' | tail -n +2) <(cat * | wc | tr -s ' '  '\n' | tail -n +2) | bc | tr -s '\n'  ' ' | sed 's/^/   /' | sed 's/$/\ /'
