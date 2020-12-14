import parse
import config
import expand
import json_ast

import copy

import os

TEST_PATH = "./tests/expansion"

if not config.config:
        config.load_config()
config.read_vars_file(os.path.join(TEST_PATH, "sample.env"))
print(config.config)

print("Using parser {} to parser tests from {}".format(config.PARSER_BINARY, TEST_PATH))

tests = os.listdir(TEST_PATH)
tests = [test for test in tests if test.endswith(".sh")]
tests.sort()

print("* Analysis tests ")

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
        is_safe = expand.safe_command(ast_object)
        
        if is_safe != expected_safe:
            print("{} command #{} expected {} got {}".format(test_name, i, expected_safe, is_safe))
            failures.add(test_name)

if len(failures) == 0:
    print("All {} tests passed".format(len(tests)))
else:
    print("{}/{} tests failed: {}".format(len(failures), len(tests), failures))

print("\n* Expansion tests")

failures = set()
for test_name in tests:
    test = os.path.join(TEST_PATH, test_name)
    ast_objects = parse.shell_file_to_ir(test)

    expected_safe = test_name.startswith("safe")
    for (i, ast_object) in enumerate(ast_objects):
        try:
            cmd = expand.expand_command(ast_object, copy.deepcopy(config.config))
            print(test_name, "expanded to", json_ast.ast_to_shell(cmd))
        except expand.EarlyError as e:
            if expected_safe:
                print("Unexpected early error in safe script", test_name)
                print("\t", e)
                failures.add(test_name)
            else:
                print("Caught early error in unsafe script", test_name)
        except Exception as e:
            print("Error:", e)
            failures.add(test_name)

if len(failures) == 0:
    print("All {} tests passed".format(len(tests)))
else:
    print("{}/{} tests failed: {}".format(len(failures), len(tests), failures))
