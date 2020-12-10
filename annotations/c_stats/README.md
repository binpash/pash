This directory contains different sets of commands; the diagram below captures most of the overlap between them, but in practice there's no subset relationship.

```
  {ubuntu { extended { POSIX {gnu coreutils} } } }
```

Each one of these sets are well-known, but a typical user setup probably contains more commands (see `ubuntu`).
To view set difference, run the following (optionally piped into `wc -l` to count):

```
../../scripts/set-diff.sh coreutils posix
```

To find all commands in an environment (e.g., to generate `ubuntu`), run `compgen -c`.

