import parse
import json
import wrapper
import config
import ast_to_ir

import os

TEST_PATH = "./tests/expansion"

print("Using parser {} to parser tests from {}".format(config.PARSER_BINARY, TEST_PATH))

tests = os.listdir(TEST_PATH)
tests.sort()

def safety(b):
    if b:
        return "safe"
    else:
        return "unsafe"

failures = set()
for test_name in tests:
    test = os.path.join(TEST_PATH, test_name)
    ast_objects = parse.shell_file_to_ir(test)

    expected_safe = test_name.startswith("safe")
    for (i, ast_object) in enumerate(ast_objects):
        is_safe = ast_to_ir.safe_command(ast_object)
        
        if is_safe != expected_safe:
            print("{} command #{} expected {} got {}".format(test_name, i, expected_safe, is_safe))
            failures.add(test_name)

if len(failures) == 0:
    print("All {} tests passed".format(len(tests)))
else:
    print("{}/{} tests failed: {}".format(len(failures), len(tests), failures))
