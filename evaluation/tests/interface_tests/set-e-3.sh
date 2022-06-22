set -e
# individual command in a multi-command pipeline
false | :
echo passed pipeline
# part of a compound list of an 'elif'
if false; then :; elif false; then :; fi
echo passed elif
# non-subshell compound command whose exit status was the result
# of a failure while -e was being ignored
{ false && : ; }
echo passed compound-brace
for i in a; do false && : ; done
echo passed compound-for
# case x in x) false && : ;; esac
# echo passed compound-case
if :; then false && : ; fi
echo passed compound-if
cont=y; while [ $cont = y ]; do cont=n; false && : ; done
echo passed compound-while
end=n; until [ $end = y ]; do end=y; false && : ; done
echo passed compound-until