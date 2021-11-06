dotFile=set.sh.tempfile
variable="value value"

# the problem is that this returns more things (we have functions that are exported in set)
set | grep variable > $dotFile
. ./$dotFile
# set
