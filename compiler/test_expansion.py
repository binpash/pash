import parse
import config
from shell_ast import expand
import json_ast

import copy

import os
import traceback

TEST_PATH = "./tests/expansion"

if not config.config:
        config.load_config()
config.read_vars_file(os.path.join(TEST_PATH, "sample.env"))
#print(config.config)

def load_ast(file):
    return json_ast.parse_json_ast_string(parse.parse_shell(test))

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
    ast_objects = load_ast(test)

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
    ast_objects = load_ast(test)

    expanded = os.path.join(TEST_PATH, test_name.replace("sh","expanded"))
    expected_safe = os.path.exists(expanded)
    for (i, ast_object) in enumerate(ast_objects):
        try:
            cmd = expand.expand_command(ast_object, copy.deepcopy(config.config))
            got = json_ast.ast_to_shell(cmd, verbose=False)

            # ??? MMG 2020-12-17 unsure about fixing the pretty printing (which may need these backslashes!)
            got = got.replace("\\'", "'")

            if not expected_safe:
                print("Unexpected success in", test_name)
                print(got)
                failures.add(test_name)
            else:
                expected = open(expanded).read()

                if got != expected:
                    print("Expected:\n\t",expected,"Got:\n\t",got)
                    failures.add(test_name)
        except (expand.EarlyError, expand.StuckExpansion,expand.Unimplemented) as e:
            if expected_safe:
                print("Found unexpected failure in", test_name)
                print("Error:", traceback.format_exc())
                failures.add(test_name)
            else:
                print("Found expected failure in", test_name)
        except Exception as e:
            print("Error:", traceback.format_exc())
            failures.add(test_name)

if len(failures) == 0:
    print("All {} tests passed".format(len(tests)))
else:
    print("{}/{} tests failed: {}".format(len(failures), len(tests), failures))
