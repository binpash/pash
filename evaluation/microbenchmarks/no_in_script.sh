N=1000 
# BUG: seq 1 $N | sort -rn
# In order to make the above work, we need to somehow get the value of N in PaSh.
# Even exporting doesn't seems to work, even though this works normally if we have the variable in the environment script.
seq 1 1000 | sort -rn
