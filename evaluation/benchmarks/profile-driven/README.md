Contains an experiment to check the profile-driven width lowering.

If we run:

```sh
time $PASH_TOP/pa.sh --parallel_pipelines --profile_driven -d 1 -w 16 --log_file pash.log for-loop.sh
```

vs

```sh
time $PASH_TOP/pa.sh --parallel_pipelines -d 1 -w 16 --log_file pash.log for-loop.sh
```

there is an obvious execution time difference.
