```
Command:    Reduction _g_:
tr            a ++ b
tr -s         (init a) ++ (tr -s (last a) ++ (head b)) ++ (tail b)
nl            a ++ map (\e -> ( (last unzip a) + (fst e)), snd e) ) b // it's simpler than it looks
join          a ++ b
cat           a ++ b
rev           b ++ a
uniq -c       a + b
uniq          (init a) ++ (uniq (last a) ++ (head b)) ++ (tail b)
uniq -d       (init a) ++ (uniq -d (last a) ++ (head b)) ++ (tail b)

join -e       
wc            ( ((fst a) + (fst b)), ((med a) + (med b)), ((lst a) + (lst b)) )
wc -l         a + b

sort          sort -m a b
sort -h       a
```
