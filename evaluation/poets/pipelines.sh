# Unix for poets: https://www.cs.upc.edu/~padro/Unixforpoets.pdf

# FIXME: Some of these are written for BSD Unix, need to re-write for GNU/Linux
# FIXME: Could we use all books from project Gutenberg?  Possibly replace with `find`
# FIXME: Features call to a script (trigrams)

# Count words
tr -sc ’[A-Z][a-z]’ ’[\012*]’ < genesis | sort | uniq -c > genesis.hist

# Merge upper and lower counts
tr ’[a-z]’ ’[A-Z]’ < genesis | tr -sc ’[A-Z]’ ’[\012*]’ | sort | uniq -c

# Count vowel sequences
tr ’a-z’ ’[A-Z]’ < genesis | tr -sc ’AEIOU’ ’[\012*]’| sort | uniq -c

# Sort 
tr -sc ’[A-Z][a-z]’ ’[\012*]’ < genesis | sort | uniq -c | sort –nr

# Bigrams (contrary to our version, this uses intermediary files)
tr -sc ’[A-Z][a-z]’ ’[\012*]’ < genesis > genesis.words
tail +2 genesis.words > genesis.nextwords
paste genesis.words genesis.nextwords |
sort | uniq -c > genesis.bigrams

# Recursive PaSh calls to trigrams (defined below)
grep ’the land of’ genesis | sh trigram | sort -nr | sed 5q
grep ’And he said’ genesis | sh trigram | sort -nr | sed 5q

# This is part of Sec. 8, not complete

# # Find words that appear 300 or more times
# tr -sc ’[A-Z][a-z]’ ’[\012*]’ < genesis | sort | uniq -c | awk ’$1 > 300 {print $2}’
# 
# # Find palindromes in genesis (nv this can be written in a single line)
# sort -u genesis.words > genesis.types
# rev < genesis.types > genesis.types.rev
# paste genesis.types genesis.types.rev | awk ’$1 == $2’

# count consonant sequences
tr ’[a-z]’ ’[A-Z]’ < genesis | tr -sc ’BCDFGHJKLMNPQRSTVWXYZ’ ’[\012*]’ | sort | uniq -c

# Sort words in Genesis by folding case.
tr -sc ’[A-Z][a-z]’ ’[\012*]’ < genesis | sort | uniq -c | sort -f

# Sort words in Genesis by rhyming order.
tr -sc ’[A-Z][a-z]’ ’[\012*]’ < genesis | sort | uniq -c | rev | sort | rev

# Count tri-grams
tr -sc ’[A-Z][a-z]’ ’[\012*]’ < genesis > genesis.words
tail +2 genesis.words > genesis.nextwords
tail +3 genesis.words > genesis.nextwords2
paste genesis.words genesis.nextwords genesis.nextwords2 |
sort | uniq -c > genesis.trigrams

# Uppercase words,
# by token
tr -sc ’[A-Z][a-z]’ ’[\012*]’ < genesis | grep -c ’ˆ[A-Z]’
# by type
tr -sc ’[A-Z][a-z]’ ’[\012*]’ < genesis | sort -u | grep -c ’ˆ[A-Z]’

# four-letter words
tr -sc ’[A-Z][a-z]’ ’[\012*]’ < genesis | grep -c ’ˆ....$’
tr -sc ’[A-Z][a-z]’ ’[\012*]’ < genesis | sort -u | grep -c ’ˆ....$’

# words with no vowels
tr -sc ’[A-Z][a-z]’ ’[\012*]’ < genesis | grep -vi ’[aeiou]’ | sort | uniq -c

# 1-syllable words
tr -sc ’[A-Z][a-z]’ ’[ 12*]’ < genesis | grep -i ’ˆ[ˆaeiou]*[aeiou][ˆaeiou]*$’ | sort | uniq -c | sed 5q

# 2-syllable words
tr -sc ’[A-Z][a-z]’ ’[\012*]’ < genesis | grep -i ’ˆ[ˆaeiou]*[aeiou][ˆaeiou]*[aeiou][ˆaeiou]*$’ | sort | uniq -c | sed 5q

# verses with 2 or more, 3 or more, exactly 2 instances of light.
grep -c ’light.*light’ genesis
grep -c ’light.*light.*light’ genesis
grep ’light.*light’ genesis | grep -vc ’light.*light.*light’

# Count morphs in genesis
spell -v genesis | sed ’s/ .*//g’ | sort | uniq -c

# Sort words by number of syllables
tr -sc ’[A-Z][a-z]’ ’[\012*]’ < genesis | sort -u > genesis.words
tr -sc ’[AEIOUaeiou\012]’ ’ ’ < genesis.words | awk ’{print NF}’ > genesis.syl
paste genesis.syl genesis.words | sort -nr | sed 5q

# Vowel sequencies that appear >=1K times
tr -sc ’[A-Z][a-z]’ ’[\012*]’ < genesis | tr -sc ’AEIOUaeiou’ ’[\012*]’ | sort | uniq -c | awk ’$1 >= 1000’

# Bigrams that appear twice
awk ’$1 == 2 {print $2, $3}’ genesis.bigrams

# Find anagrams
rev < genesis.types > genesis.types.rev
sort genesis.types genesis.types.rev |
uniq -c |
awk ’$1 >= 2 {print $2}’

# Compare Exodus and Genesis
tr -sc ’[A-Z][a-z]’ ’[\012*]’ < genesis | sort -u > genesis.types
tr -sc ’[A-Z][a-z]’ ’[\012*]’ < exodus | sort -u > exodus.types
sort genesis.types exodus.types exodus.types | uniq -c | head

# Fix awk implementation of uniq -c (nv: unclear what this is)
awk ’$0 == prev { c++ }
$0 != prev { print c, prev
c=1
prev=$0 }
END {print c, prev}’
