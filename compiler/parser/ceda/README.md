# JSON <-> shell script (without OCaml)

## Pre-requisites

* json-c v0.15:

```
cd /pash
wget https://s3.amazonaws.com/json-c_releases/releases/json-c-0.15.tar.gz
tar zxf json-c-0.15.tar.gz
cd json-c-0.15
mkdir build
cd build/
cmake -DCMAKE_INSTALL_PREFIX=/pash/json-c-0.15/install ../
make install
```

## JSON -> shell script

### Usage

It has the same usage as `json_to_shell` i.e.,
`./json_to_shell2 JSON_FILE_NAME` or `./json_to_shell2 < JSON_FILE_NAME`

e.g.,

```
$ cat /pash/evaluation/hello-world.sh | ../parse_to_json.native | ./json_to_shell2
if [ $(uname) = "Darwin" ]; then a=/usr/share/dict/web2; else a=/usr/share/dict/words; fi
if [ -f ${a} ]; then cat ${a} ${a} ${a} ${a} ${a} ${a} ${a} ${a} | grep "\\(.\\).*\\1\\(.\\).*\\2\\(.\\).*\\3\\(.\\).*\\4" | wc -l; else echo "Dictionary file ${a} not found.."; fi
```

### Testing

#### Testing a single script

```
sh test_JSON_to_shell2.sh SOME_SCRIPT_FILE
```

This applies the OCaml implementation of `parse_to_json` to the specified script file, then compares the OCaml `json_to_shell` against this re-implementation.

Output:
* `REF_ABORT_1` or `REF_ABORT_2` means the shell script cannot be parsed to or from JSON by the reference implementation. 
* `ABORT` or `FAIL` means there's a bug in Thurston's code.
* `PASS` is good.

#### Testing all the scripts in /pash/

This includes `test_JSON_to_shell2.sh`!

```
make testsA-all
```

All shell scripts that the OCaml implementation works on are regenerated, byte-for-byte identical:
```
    234 PASS
     33 REF_ABORT_1
```

### Known Bugs

* `fresh_marker` for heredocs. This is really obscure and a pain to implement in C. For real-world, non-adversarial settings, just change the marker from "EOF" to some random text.
* Not Python

## WORK-IN-PROGRESS: Shell script -> JSON

### Usage

It has the same usage as `parse_to_json` i.e.,
`./parse_to_json2 SCRIPT_FILE_NAME` or `./parse_to_json2 < SCRIPT_FILE_NAME`

e.g.,
```
$ ./parse_to_json2 /pash/evaluation/hello-world.sh
```

### Testing

#### Testing a single script

```
sh test_parse_to_JSON2.sh SOME_SCRIPT_FILE
```

This applies the OCaml implementation of `parse_to_json` to the specified script file, and compares the output against this re-implementation.

Output:
* `REF_ABORT_1` means the shell script cannot be parsed to JSON by the reference implementation. 
* `ABORT` or `FAIL` means there's a bug in Thurston's code.
* `PASS` is good.

#### Testing all the scripts in /pash/

```
make testsB-all
```

Current progress (come back soon for updates!):
```
     49 ABORT
     29 FAIL
    174 PASS
     33 REF_ABORT_1
```
