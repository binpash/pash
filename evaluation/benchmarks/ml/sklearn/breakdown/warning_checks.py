import sys
from pickle import dumps, load

reg = load(sys.stdin)
C_, penalty = reg.warning_checks()

sys.stdout(dumps(C_))
sys.stdout(dumps(penalty))