#!/bin/bash

## Initialize the necessary temporary files
file1=$(mktemp)
file2=$(mktemp)
file3=$(mktemp)
file4=$(mktemp)
file5=$(mktemp)
file6=$(mktemp)
file7=$(mktemp)

export LC_ALL=C

# Commit history in the form of ascending Unix timestamps, emails
git log --pretty=tformat:'%at %ae' | 
awk 'NF == 2 && $1 > 100000 && $1 < '`date +%s` | 
sort -n > "$file1"

# Calculate number of committers
awk '{print $2}' "$file1" | 
sort -u | 
wc -l > "$file2"
cp "$file2" "$file3"
cp "$file2" "$file4"

# Calculate last commit timestamp in seconds
tail -1 "$file1" | 
awk '{print $1}' > "$file5"

# Calculate first commit timestamp in seconds
head -1 "$file1" | 
awk '{print $1}' >> "$file5"

# Gather last and first commit timestamp and compute the difference in days
cat "$file5" | 
tr '\n' ' ' | 
awk '{print int(($1 - $2) / 60 / 60 / 24)}' > "$file5"

sort -k2 "$file1" > "$file6"

# Place committers left/right of the median according to the number of their commits
awk '{print $2}' "$file1" | 
sort | 
uniq -c | 
sort -n | 
awk -v committers1="$file2" '
BEGIN {
    while ((getline NCOMMITTERS < committers1) > 0) {}
    l = 0; r = NCOMMITTERS;
}
{print NR % 2 ? l++ : --r, $2}' |
sort -k2 > "$file7"

# Join committer positions with commit timestamps based on committer email
join -j 2 "$file6" "$file7" | 
sort -k 2n > "$file6"

# Create portable bitmap
{
    echo 'P1'
    {
        cat "$file3"
        cat "$file5"
    } | 
    tr '\n' ' ' | 
    awk '{print $1, $2}'
    
    perl -na -e '
    BEGIN {
        open(my $ncf, "<", "'"$file4"'");
        $ncommitters = <$ncf>;
        @empty[$ncommitters - 1] = 0; @committers = @empty;
    }
    sub out {
        print join("", map($_ ? "1" : "0", @committers)), "\n";
    }

    $day = int($F[1] / 60 / 60 / 24);
    $pday = $day if (!defined($pday));

    while ($day != $pday) {
        out();
        @committers = @empty;
        $pday++;
    }

    $committers[$F[2]] = 1;

    END { out(); }
    ' "$file6"
} | 
pgmmorphconv -erode <(
cat <<EOF
P1
7 7
1 1 1 0 1 1 1
1 1 0 0 0 1 1
1 0 0 0 0 0 1
0 0 0 0 0 0 0
1 0 0 0 0 0 1
1 1 0 0 0 1 1
1 1 1 0 1 1 1
EOF
) | 
tee | 
{
    # Full-scale image
    pnmtopng >large.png
    # A smaller image
    pamscale -width 640 | 
    pnmtopng >small.png
}
