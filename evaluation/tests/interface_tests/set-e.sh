set -e
# part of a compound list of a 'while', 'until' or 'if'
while false; do break; done
echo passed while
until false; do break; done
echo passed until
if false; then :; fi
echo passed if
# any command of an AND-OR list other than the last
false && :
echo passed AND list
false || :
echo passed OR list
: && false || false && :
echo passed AND-OR list
# part of a pipeline preceded by the '!' reserved word
! false
echo passed negated pipeline