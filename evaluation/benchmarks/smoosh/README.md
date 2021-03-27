## Running Smoosh Tests

1. Clone the [Smoosh repository](https://github.com/mgree/smoosh.git)
2. Apply the patch `smoosh_problematic_tests.patch` which comments out some tests that lead to space issues
3. Make sure to export `$PASH_TOP` and the other final lines in `scripts/install.sh`
4. When in the smoosh top directory, run:

```
TEST_TIMEOUT=10 TEST_DEBUG=1 TEST_SHELL="${PASH_TOP?No PASH_TOP}/pa.sh --width 2 --log_file /tmp/log_file" make -C tests veryclean
TEST_TIMEOUT=10 TEST_DEBUG=1 TEST_SHELL="${PASH_TOP?No PASH_TOP}/pa.sh --width 2 --log_file /tmp/log_file" make -C tests
```

Running `veryclean` before the tests is important since some things are cached.

If you just want to see a summary run it without `TEST_DEBUG=1`.
