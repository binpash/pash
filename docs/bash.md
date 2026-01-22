# Bash Mode in PaSh

## Usage

Add `--bash` to the pash call. For example:
```sh
pash --bash ./evaluation/intro/hello-world-bash.sh
```

Bash mode is experimental -- bug reports are welcome!

## Design

PaSh originally targeted POSIX-compliant scripts through libdash and shasta.
The `--bash` flag enables the alternative libbash front-end.

It works similarly: libbash uses bash as a library to parse the script,
which shasta interprets as an in-python AST for pash to analyze.
Shasta unifies both dash and bash ASTs, with additional bash-only nodes and fields.
However, shell expansion for `--bash` is implemented through a bash mirror server, instead of in-python like in the libdash mode.

Pash doesn't support any bash-specific features, but it could be extended to in the future.

## Limitations

Libbash does not parse shell metacharacters into ArgChars. Example:
```sh
echo "$(echo)"
```

Libdash, through shasta, parses this as:
```
Command:
- echo
- QArgChar
   - BArgChar
       - Command: echo
```
Libbash, however, does not have that granularity:

```
Command:
- echo
- "$(echo)"
```

This makes proper shell expansion -- necessary for parallelization -- impossible in python,
so, potentially expandable words are echoed to a bash mirror process to expand.
However, the bash server has to be conservative. Notable limitations:
- Arguments containing `*` leads to the command becoming unparallizable, since it could be a glob.
- For loops can't be parallelized because the looping variable would look like a necessary expansion.

A future project is to add ArgChar parsing and support to libbash, shasta and sh-expand.

