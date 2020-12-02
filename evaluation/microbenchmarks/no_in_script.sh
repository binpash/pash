N=100
# BUG: seq 1 $N | sort -rn
# In order to make the above work, we need to somehow get the value of N visible in the PaSh runtime python process.
seq 1 100 | sort -rn
