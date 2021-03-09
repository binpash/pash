#!/bin/sh


input_script='/pash/compiler/parser/libdash/ltmain.sh'


if [ $# -eq 1 ]
then
    input_script="$1"
fi


echo "Input script: $input_script"
echo

echo "OCaml (dash C AST -> libdash OCaml AST -> JSON -> Pash Python AST -> JSON -> shell:"
time (../parse_to_json.native "$input_script" > /tmp/json.$$; cat /tmp/json.$$ | ../json_to_shell.native | md5sum)
echo

echo "C (dash C AST -> libdash C AST -> JSON -> Pash Python AST -> JSON -> shell):"
time (./parse_to_json2 "$input_script" > /tmp/json.$$ 2>/dev/null; cat /tmp/json.$$ | ./json_to_shell2 | md5sum)
echo

echo "Python (dash C AST -> libdash C AST -> JSON -> Pash Python AST -> JSON -> shell):"
time (python3 ./parse_to_json2.py "$input_script" > /tmp/json.$$ 2>/dev/null; cat /tmp/json.$$ | python3 ./json_to_shell2.py | md5sum)
echo

echo "Python (dash C AST -> Pash Python AST -> shell):"
time (python3 ceda_rt.py "$input_script" | md5sum)
echo

