# PaSh Issue #741: Quoted Commands Are Not Recognized as Parallelizable

## Status

This document describes the solution and regression coverage for
[binpash/pash#741](https://github.com/binpash/pash/issues/741).

## Issue summary

PaSh recognizes commands through its annotation library before constructing and
optimizing a dataflow graph. A command written without quotes was recognized:

```bash
cat README.md | tr "A-Z" "a-z"
```

The shell-equivalent form with quoted executable names was not recognized:

```bash
"cat" README.md | "tr" "A-Z" "a-z"
```

Before this fix, running the quoted form with
`--assert_all_regions_parallelizable` failed because PaSh reported that zero of
one regions were parallelized. The unquoted form succeeded.

## Root cause

Shell quotes control word expansion; they are not part of the resulting
executable name. Shasta's AST preserves quoted portions as `QArgChar` nodes.

During annotation lookup, PaSh formatted the entire command-name AST with the
general-purpose `format_arg_chars` function. `QArgChar.format()` intentionally
reintroduces double quotes for shell rendering, so the annotation lookup received
`"cat"` rather than `cat`. Since annotations exist for `cat`, not for a command
literally named `"cat"`, the command was treated as unannotated and therefore
side-effectful.

## Solution

The fix introduces `remove_quotes_expanded_arg_chars`, a command-name-specific
formatter that:

1. Accepts normal (`CArgChar`) and escaped (`EArgChar`) characters.
2. Recursively traverses quoted (`QArgChar`) nodes without emitting their
   syntactic quote wrappers.
3. Rejects any other node type with `ValueError`, because command names must be
   fully expanded before annotation lookup.
4. Is used only to obtain the executable name for annotation lookup. Arguments,
   operands, flags, and emitted shell code continue using the existing formatter.

Restricting normalization to the executable name is important: removing quotes
from operands could change shell semantics through field splitting, globbing, or
empty-argument handling.

## Changed files

- `src/pash/compiler/util.py`
  - Adds the recursive `remove_quotes_expanded_arg_chars` helper.
- `src/pash/compiler/annotations_utils/util_parsing.py`
  - Uses the helper when looking up annotations for a command name.
- `evaluation/tests/quoted_cmd.sh`
  - Adds fully quoted, single-quoted, and mixed-quote command-name regressions.
- `evaluation/tests/test_evaluation_scripts.sh`
  - Registers `quoted_cmd` as a pipeline microbenchmark that must be
    parallelizable.

## Behavioral coverage

The local validation covers:

- Unquoted command names: `cat`
- Fully double-quoted command names: `"cat"`
- Fully single-quoted command names: `'cat'`
- Adjacent mixed segments: `c"a"t`
- Escaped characters: `c\at`
- Command names produced by variable expansion: `cmd=cat; "$cmd" ...`
- Quoted operands, ensuring command-name normalization does not alter arguments
- Nested quoted AST nodes
- Empty quoted AST nodes
- Rejection of unexpanded variable AST nodes
- Rejection of unknown quoted executables
- Rejection of executable names containing literal quote characters

## Reproduction

From the repository root with the local Python environment active:

```bash
export PASH_TOP="$PWD/src/pash"
export PATH="/opt/homebrew/bin:$PWD/.venv/bin:$PATH"

# Expected to succeed both before and after the fix.
pash --dry_run_compiler --assert_all_regions_parallelizable -d 1 \
  -c 'cat README.md | tr "A-Z" "a-z"'

# Failed before the fix and succeeds after it.
pash --dry_run_compiler --assert_all_regions_parallelizable -d 1 \
  -c '"cat" README.md | "tr" "A-Z" "a-z"'
```

## Regression test

The focused regression can be executed without parallel runtime optimization:

```bash
export PASH_TOP="$PWD/src/pash"
export PATH="/opt/homebrew/bin:$PWD/.venv/bin:$PATH"

/opt/homebrew/bin/bash evaluation/tests/quoted_cmd.sh > /tmp/quoted-sequential.out
pash --no_optimize --assert_all_regions_parallelizable -d 1 \
  evaluation/tests/quoted_cmd.sh > /tmp/quoted-pash.out
cmp /tmp/quoted-sequential.out /tmp/quoted-pash.out
```

A successful run exits with status zero and produces byte-for-byte identical
output.

## Test platforms

The checkout is hosted on Apple Silicon macOS. PaSh's optimized runtime cannot
currently be built directly there because the upstream `dgsh-macho.s` fails to
assemble with the installed Apple toolchain. Host-side validation therefore uses
`--dry_run_compiler` for optimized compilation and `--no_optimize` for semantic
output comparison.

Optimized end-to-end validation was additionally performed in a clean, native
ARM64 Ubuntu 24.04 container using Docker with Colima. All runtime helpers,
including `dgsh-tee`, were rebuilt inside Linux and verified as ARM64 ELF
executables before tests were run. The repository was mounted read-only while
the test tree was copied into the container, so the container could not modify
the working checkout.

## Local validation results

The complete issue-specific validation matrix was executed twice from a clean
test-output directory. Both runs completed with `VALIDATION_FAILURES=0` and exit
status zero. Each run included:

- Seven successful positive end-to-end compiler cases
- Three correctly rejected negative compiler cases
- Byte-for-byte Bash-versus-PaSh output comparison
- Seven direct AST normalization cases
- Rejection of an unexpanded variable AST node
- Quoted-versus-unquoted annotation-parser equivalence
- Bash syntax validation
- Python byte-compilation
- Black formatting validation for both modified Python modules
- `git diff --check`

The new `quoted_cmd.sh` regression was then run in Linux in every configuration
used by the compiler test harness. All four runs exited successfully, satisfied
`--assert_all_regions_parallelizable`, and produced output byte-for-byte
identical to Bash:

- `--width 2 --bash`
- `--width 8 --bash`
- `--width 2` (Dash backend)
- `--width 8` (Dash backend)

The complete Linux issue-specific matrix (positive cases, expected rejections,
all four semantic configurations, AST checks, and annotation-parser checks) was
run twice with `LINUX_VALIDATION_FAILURES=0`. It was run once more after the
final mechanical Black formatting change, again with zero failures.

The repository-wide `scripts/run_tests.sh` was also invoked in the clean Linux
environment. Its optimized intro suite passed 3/3. The existing interface suite
reported 23/43 passing, with failures across shell arguments, redirections,
`set -e`, traps, globbing, and IFS behavior rather than the quoted-command path.
The compiler suite then stalled in the existing `comm-par-test` and
`comm-par-test2` cases (`comm -23` remained blocked). Terminating their orphaned
processes also terminated the aggregate harness, so no complete repository-wide
pass is claimed. The issue-specific test was run independently afterward with
per-run timeouts and passed all four configurations as described above.

A second repository-wide run used `--no_optimize` to bypass the missing parallel
runtime. Its completed intro cases passed, but the existing `hello-world-bash`
case did not terminate on this macOS host, so that run was also stopped rather
than reported as a complete suite pass.

Finally, source-distribution generation completed, but wheel generation reached
the project's runtime build step and failed because `eager_lib.c` references the
Linux-oriented `std_copy` and `send_file` helpers that are unavailable with the
Apple toolchain. No packaging error involved the files changed for issue #741.

## Review checklist before any commit

- [ ] Review the implementation and regression test.
- [ ] Confirm all local validation results.
- [ ] Investigate or confirm the unrelated Linux interface failures and hanging
      `comm` tests in PaSh's CI environment if a repository-wide green run is
      required.
- [ ] Decide whether this documentation file should be included in the eventual
      commit or retained only as local review notes.
- [ ] Commit only after explicit approval.
- [ ] Push or open a pull request only after separate explicit approval.
