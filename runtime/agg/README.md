Currently aggregators are WIP. The new ones are in `cpp/bin`. They are automatically built during `setup_pash.sh` and the unit tests in `cpp/tests` are run during `run_tests.sh`. The interface is like the following:

```sh
aggregator inputFile1 inputFile2 outputFile args
```

Where `args` are the arguments that were passed to the command that produced the input files.
