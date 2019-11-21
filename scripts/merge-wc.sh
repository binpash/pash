paste -d '+' <(cat * | wc | tr -s ' '  '\n' | tail -n +2) <(cat * | wc | tr -s ' '  '\n' | tail -n +2) | bc
