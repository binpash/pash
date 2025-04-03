This directory contains different groups of commands; the diagram below captures most of the overlap between them, but in practice there's no real subset relationship.

```
  {ubuntu { extended ( POSIX {gnu )coreutils} } } }
```
Each one of these sets are well-known The focus is on POSIX and GNU Coreutils. The study and annotations of other commands is welcome!

To view set difference between classes of commands, run the following script from `pash`'s benchmarks (and count them by piping to `wc -l`):

```
../../scripts/set-diff.sh coreutils posix
```

To find all commands in an environment (e.g., to generate `ubuntu`), run `compgen -c`.

