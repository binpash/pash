# JSON <-> shell script (without OCaml)

## JSON -> shell script

### Pre-requisites

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

### Usage

It has the same usage as `json_to_shell` e.g.,

```
$ ./json_to_shell2 hello_world.json
if [ $(uname) = "Darwin" ]; then a=/usr/share/dict/web2; else a=/usr/share/dict/words; fi
if [ -f ${a} ]; then cat ${a} ${a} ${a} ${a} ${a} ${a} ${a} ${a} | grep "\\(.\\).*\\1\\(.\\).*\\2\\(.\\).*\\3\\(.\\).*\\4" | wc -l; else echo "Dictionary file ${a} not found.."; fi
```

It can also read from standard input instead.

### Testing

#### Testing a single script

```
sh test_JSON_to_shell2.sh SOME_SCRIPT_FILE
```

This applies the OCaml implementation of `parse_to_json` to the specified script file, then compares the OCaml `json_to_shell` against this re-implementation.

Note that some shell scripts cannot be parsed to or from JSON by the reference implementation; these are denoted by the errors `REF_ABORT_1` or `REF_ABORT_2`. `ABORT` or `FAIL` means there's a bug in Thurston's code. `PASS` is good.

#### Testing all the scripts in /pash/
```
make tests-all
```

All shell scripts that the OCaml implementation works on are regenerated, byte-for-byte identical:
```
    234 PASS
     33 REF_ABORT_1
```

### Known Bugs

* `fresh_marker` for heredocs. This is really obscure and a pain to implement in C. For real-world, non-adversarial settings, just change the market from "EOF" to some random text.
* Not Python

## Shell script -> JSON (TODO)
