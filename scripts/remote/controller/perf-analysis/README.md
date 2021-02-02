This subsystem analyzes performance test results, and summarizes them
for consumption by either humans or other programs.

Each performance test computes timing information, but tests
themselves have variants defined by the PaSh paper (sequential or
distributed in some way).

Information about a test impacts the name, location, and content of
that test's results file.  To summarize and compare findings from such
files, you need to know the test's name, variant, and DFG width. If
the test has additional options (namely input size), you may also have
to know the subdirectory in which the file appears.


# Building a Report
[reporting]: #building-a-report

Assuming a populated `results` directory, `reportPerfSuite` from the
`report.js` module will show the most helpful timing information it
can in terms of the arguments. The output is meant for human
consumption. If you want to deal with the raw data from the file, then
see [here][reading].

```javascript
> const { reportPerfSuite } = require('./report.js');
> console.log(reportPerfSuite('results', ['diff', 'wf', 'sort'], 2, 'distr'))

╔════════════════╤════════════════╤════════════════╗
║ diff           │ sort           │ wf             ║
╟────────────────┼────────────────┼────────────────╢
║ 959.04s, x2.64 │ 717.36s, x1.65 │ 742.16s, x1.78 ║
╚════════════════╧════════════════╧════════════════╝
```

PaSh execution time for some shell program is shown in seconds. The
`x`-prefixed scalar is a speedup value that indicates how fast PaSh
was going compared to the (sequential) shell.

If no `*_seq.time` files are present to add that context, then the
speedup value is not shown.


# Reading Performance Test Results
[reading]: #reading-performance-test-results

`parse.js` reads the performance test result files and converts them
to an in-memory representation. This is useful for things like CI,
which needs to track metrics over different revisions.

```javascript
const path = require('./path');
const p = require('./parse.js');
const resultsDir = '/path/to/results';

// Creates metadata for a single file
const fileName = p.makePerfFileName('wf', 2, 'distr');
const data = p.analyzePerfFile(path.join([resultsDir, fileName]));

// A summary is just a human readable string that shows
// the most timing information it can.
console.log('Summary: ', p.summarizePerfData(data));

// Creates metadata for a set of files that covers the given tests, for
// a fixed DFG width and test variant.
console.log(p.analyzePerfSuite(resultsDir, ['wf', 'diff', 'sort'], 2, 'distr'));
```

*NOTE:* `parse.js` also exports `fileContentCache`, which holds
lazily-loaded perf file contents in memory for a limited time. Node
will wait on these timeouts at the end of its process.  To avoid
waiting at the end of the process, clear the cache with
`fileContentCache.clear()`.
