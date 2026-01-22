# Runtime Support
Quick Jump: [Stream Splitting](#stream-splitting) | [Eager Stream Polling](#eager-stream-polling) | [Cleanup Logic](#cleanup-logic) | [Aggregators](#aggregators)

PaSh includes a small library of runtime primitives supporting the runtime execution of parallel scripts emitted by the compiler.

## Stream Splitting

The PaSh compiler inserts `split` nodes to expose parallelism when parallelizable nodes only have one input.

## Eager Stream Polling

To overcome the laziness challenges outlined in Sec. 5, PaSh inserts and instantiates `eager` nodes on streams.

## Cleanup Logic

PaSh contains cleanup logic for dealing with dangling FIFOs.
This is implemented in `wait_for_output_and_sigpipe_rest.sh`.

## Aggregators

There is a small custom aggregator library provided in [agg/py/](agg/py/).
These aggregators are used to merge partial results from the parallel scripts.

For example, the aggregator `wc.py` can merge results from partial `wc` running in parallel.
To confirm what the aggregator does, call it as follows.

```shell 
$PASH_TOP/runtime/agg/py/wc.py <(echo -e '1\n2\n3' | wc -l) <(echo -e '1\n2\n3' | wc -l)
```

Internally, the aggregator looks roughly like this:

```Python
#!/usr/bin/python
import sys, os, functools, utils

def parseLine(s):
  return map(int, s.split())

def emitLine(t):
  return [" ".join(map(lambda e: str(e).rjust(utils.PAD_LEN, ' '), t))]

def combiner(a, b):
  if not a:
    return b
  az = parseLine(a[0])
  bz = parseLine(b[0])
  return emitLine([ (i+j) for (i,j) in zip(az, bz) ])

utils.help()
res = functools.reduce(agg, utils.read_all(), [])
utils.out("".join(res))
```

The core of the aggregator, function `combiner`, is binary (i.e., takes two input streams).
The `reduce` function lifts the aggregator to arity _n_, i.e., the number of the incoming edges—each of which feeds the aggregator with the results of running a parallel `wc` instance.
This allows developers to think of aggregators in terms of two inputs, but generalize their aggregators to many inputs.
Utility functions such as `read` and `help`, common across our aggregator library, deal with error handling when reading multiple file descriptors, and offer a invocation flag `-h` that demonstrates the use of each aggregator.

PaSh’s library currently several aggregators, many of which are usable by more than one command or flag. For example, the aggregator shown above is shared among `wc`, `wc -lw`, `wc -lm` etc.


