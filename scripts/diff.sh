# https://crashingdaily.wordpress.com/2008/03/06/diff-two-stdout-streams/
diff -B <( sort A | tr [:lower:] [:upper:] ) <( sort B | tr [:lower:] [:upper:])
