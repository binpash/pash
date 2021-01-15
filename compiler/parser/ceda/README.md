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
make -j 8 install
```

### Usage

It has the same usage as `json_to_shell` e.g.,

```
./json_to_shell2 hello_world.json
if [ $(uname) = "Darwin" ]; then a=/usr/share/dict/web2; else a=/usr/share/dict/words; fi
if [ -f ${a} ]; then cat ${a} ${a} ${a} ${a} ${a} ${a} ${a} ${a} | grep "\\(.\\).*\\1\\(.\\).*\\2\\(.\\).*\\3\\(.\\).*\\4" | wc -l; else echo "Dictionary file ${a} not found.."; fi
```

* Testing:

```
make tests-all
```

This applies the OCaml implementation of `parse_to_json` to all the shell scripts in `/pash/`, then compares the OCaml `json_to_shell` against this re-implementation.

Note that some shell scripts cannot be parsed to JSON (denoted by the error `REF_ABORT_1`); all other shell scripts are regenerated, byte-for-byte identical:

Expected output:
```
    234 PASS
     33 REF_ABORT_1
```

### Known Limitations

* `fresh_marker` for heredocs. This is really obscure and a pain to implement in C. For real-world, non-adversarial settings, just change the market from "EOF" to some random text.

## Shell script -> JSON (TODO)
