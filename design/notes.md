

## Meeting 

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

Thinking about the next version (and expansion), we might need to hook directly into an existing shell.
My hope is that we can avoid being shell-specific and remain shell agnostic.
We have discussed two obvious ways:

A. Expose an interface where we call into the original shell (responsible for the global).
The shell expands a bit and then gives us the result, and maintains all the state.
Our compiler can call any shell we want, after all---including the user's original shell.

B. Use commands that invert `pash`'s evaluation order, calling `pash` after expressions have been expanded.
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

Probably the technique used at the end will require combining these solutions, and possibly invent a few more:-)



