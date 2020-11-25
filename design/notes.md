

## Meeting Notes Oct. 29 2020

There are several different modes: char mode, word mode, line, bytes
`IFS` might be important, as is `cat -b` etc., `gunzip`, `pad` (parser combinator?)

Breaking lazy semantics is dangerous, esp. for utilities such as `yes`; as such `eager` needs to be bounded and block after a specific amount
(nv: my intuition is that there's an impossibility theorem somewhere here---you can't predict how much to buffer, there will always be cases where any arbitrary limit will fail.)

There's also an implicit cost model here not made explicit (nv: not quite, but there should have been as plots in Fig.10 show).
It could be like: enumerate every DFG up to a width, let me give a cost to this.
We could also do it in one sort --- this is where the JIT setting might be helpful

`stat` is cheating; but `curl` and `wget` will tell you how many bytes the downloaded, could this be used for further optimizations?

In terms of experiments, it would be interesting to see whether these are replica table across environments and setups (also: maybe we need a better eval suite?).
See, for example, a post talking about the fact that essentially VM warmup does not exist (it's not a phenomenon); something similar could be happening here.
Put differently, to what extend are these optimizations pessimizastions.

Maybe we need an eval suite of "portable" shell scripts (but also maybe not on the critical path).

### Design Alternatives

Thinking about the next version (and expansion), we might need to hook directly into an existing shell.
It would be interesting (but not necessary) to remain shell agnostic, i.e. not extending an existing shell.
Two possible alternatives follow. Probably the technique used at the end will require combining these solutions, and possibly invent a few more:-)

#### Tight integration with a shell

Expose an interface where we call into the original shell (responsible for the global).
The shell expands a bit and then gives us the result, and maintains all the state.
Our compiler can call any shell we want, after all---including the user's original shell.

#### Preprocessing AOT pass that prepends PaSh

Use commands that invert `pash`'s evaluation order, calling `pash` after expressions have been expanded.
This includes an AOT components that prefixes every command _κ_ with, say, `pash-call <ID>`_κ_, which achieves two things.
First, it allows the shell to expand arguments before calling `pash-call` (including, for example, cases where the command itself is the result of an expansion).
Second, by generating IDs, it allows `pash` to understand where the current execution fragment (during JIT execution) is within the larger execution tree.
Here is an example:

```
echo $ONE | wc $FLAGS

# first pash AOT pass =>

pash-call ID1234 echo $ONE | pash-call ID4245 wc $FLAGS
```

Or maybe:

```
pash-pipe PIPE0
pash-call ID1234 - PIPE0 -- echo $ONE
pash-call ID2345 PIPE0 - -- wc $FLAGS
pash-run ID2345
```

Since we'll want to be able to, e.g., split `PIPE0`.

Possible issue: A possible scalability issue with this approach is that _if_ we want to perform other kinds of optimizations (not just the parallelization ones), e.g. replacing `;` with `&` if there are no dependencies (see Emery's work), we would need to replace every construct in the script (`;`, `&&`, `|`, `&`, ...) with a call to PaSh.

Question: Is there any way that the necessary structural information (pipes `|`, sequence operators `;`, ...) may be hidden AOT (leading to an "incomplete" transformation)?

Question: Do we actually need to separate each construct in a different call to PaSh or can we just call PaSh once with the complete AST of the target subprogram (as we do now) only forcing string expansion before the call?

Example (that might not be correct):

```sh
ls "${DIR}/*" | wc $FLAGS

##     |
##     v

temp1="${DIR}/*"
temp2=$FLAGS
pash-compile ID1 $temp1 $temp2
```

It depends on the string expansion, though. Consider the following
(hare-brained) example:

```sh
cat notes.md | { count=0; while read line; do : $((count += 1)); done; echo $count ; }
```

It's not possible to expand `$count` correctly early. It can be even
worse:

```sh
cat notes.md | 
  { count=0
    while read line
    do 
      : $((count += 1))
      printf "%05d\t%s\n" "$count" "$line"
    done
  }
```

Neither `$count` nor `$line` can be expanded early. On the plus side,
these variables are totally local to this mediocre implementation of
`cat` (or, as above, `wc -l`). Since only builtins can affect the
environment, we should be able to determine in most situations which
variables can be expanded and which can't.

This analysis isn't just about parameter expansion---globbing gets hit too.

```sh
cmd1 | { read x; cd $x; echo *; }
```

### Issue with shell state management

An issue with all proposed solutions (and our existing prototype) is that any changes to the shell state that are performed in a dataflow graph (a parallel script that PaSh generates and runs), are not visible in the original shell that would have executed this code. Such changes include changes to values of environment variables.

To address this we would need to either execute the parallel (optimized) script in the original shell, or manage to copy the state changes to the original shell after execution.

Some of these issues are mitigated by subshell spawning. Asynchronous
jobs and pipelines are run in subshells, which are forked from the
main shell and can't affect its state. There's one subtle exception:
POSIX shells may or may not fork a subshell for the last job of a
synchronous pipeline. That is, it is POSIX compliant for the following
program to set `x` to 5 or 10:

```sh
x=10; echo 5 | read x; echo $x
```

Bash, dash, and yash will set it to 10, i.e., the `read` occurs in a
subshell. Zsh sets it to 5.

### Safe early expansions

It's "safe" to expand something early when it either (a) won't affect
the shell state, or (b) will stop execution. (Assuming `set -o
pipefail` if you're in a pipeline.)

 0. NONPOSIX: Brace expansion
    * `{w1,w2,...}` is safe if `wi` are safe
 1. Tilde
    * always safe
 2. Parameter
    * `$x` and `${#x}` are safe
    * `${x?w}` is safe
    * `${x-w}`, `${x+w}` `${x%w}`, (and `%%`, and `#`, and `##`) are safe if `w` is safe
    * `${x=w}` is safe when `x` is already assigned, otherwise UNSAFE
    * NONPOSIX: `${x:off:len}` is safe when `off` and `len` are safe
    * NONPOSIX: `${!w}` is safe when `w` is safe
    * NONPOSIX: `${x/pat/str/}` is safe when `pat` and `str` are safe
 3. Arithmetic
    * operations are safe
    * `+=` and `=` and family are UNSAFE
    * NONPOSIX: `++` and `--` are UNSAFE
    * `op="+=1"; $((x $op))` is UNSAFE
 4. Command
    * `$(w)` is safe if `w` is safe
 5. Field splitting
    * always safe
 6. Pathname
    * safe so long as we haven't run `cd` in our block that we can't yet account for, e.g.
      `foo | { read x; cd $x; echo *; }` is UNSAFE
 7. Quote removal
    * always safe
    
Approach:

  Write a safety analysis that walks over the `arg list`
  (cf. `compile_arg_char` from `compiler/ast_to_ir.py`).
  
  A `safe_expand` procedure walks over things and expands things that
  it deems safe.
