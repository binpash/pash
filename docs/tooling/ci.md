
## Continuous Integration

A text-based continuous integration (CI) sever is set up at [pash.ndr.md](http://pash.ndr.md).
The focus of CI is to monitor correctness, not performance: 
  it runs all of PaSh's benchmarks, with small inputs, and compares results with the sequential execution.
It additionally runs and reports on the tests from the Smoosh suite.

#### Summary 

To get the summary of the latest 5 builds:

```sh
curl -s pash.ndr.md | head
```

This will output a listing of the following form; the first row does not exist in the output, and is used only to show what the output means:

```
Date       Time      Commit   Commit Message                        Tests: PaSh | Smoosh  Time
2020-12-11;12:20:18  0277819  First version CI                          144/144 | 92/174  0s
2020-12-12;09:02:13  7670562  Merge branch 'master' of github.com:a...  144/144 | 78/174  0s
2020-12-12;09:42:59  c6521b9  Update webhook                            144/144 | 78/174  0s
2020-12-12;10:33:33  f821b4b  Merge branch 'master' of github.com:a...  144/144 | 77/174  0s
```

As this is an append-only log, it may feature copies of jobs run for the same commit --- for example, when testing was run twice for the same commit.

#### Single-build Info

Getting more information about a build is achievable via:

```sh
curl -s pash.ndr.md/f821b4b
```

This will output the full log up until the current time. That is, visiting the same URL after a few seconds (while the specific job is running) will report more output.

When the job is complete, this endpoint will permanently host the output for the job. If the job is re-run, the endpoint will host only the latest log---i.e., the old log is not kept.

## Job/Status Endpoints

To launch jobs, one can run

```sh
curl -s ctrl.pash.ndr.md/<job>
```

This will launch the job, if there is no other job running; the full test suite (see job `ci` below) takes about 10 minutes, and thus there is rarely an existing job running when a new commit is pushed to the repo.
Possible `job` endpoints include:

* `ci`: Runs the tests and produces the reports mentioned above
* `pkg` Packages the latest commit as an archive available at `up.pash.ndr.md`.
* `now` Reports on any job executing currently

The `ci` job is run on every push event from GitHub, and `pkg` is run as as part of `ci`.


