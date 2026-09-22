import os
import sys

import pytest
from shasta.ast_node import (
    AArgChar,
    BArgChar,
    CArgChar,
    EArgChar,
    QArgChar,
    VArgChar,
)

sys.path.insert(
    0,
    os.path.join(os.path.dirname(__file__), "..", "src", "pash", "compiler"),
)

from util import format_arg_chars, remove_quotes_expanded_arg_chars


def literal(string):
    """Build the arg_char list for an unquoted literal word."""
    return [CArgChar(ord(char)) for char in string]


def quoted(arg_chars):
    """Wrap an arg_char list in a single quoted group."""
    return [QArgChar(arg_chars)]


def test_unquoted_name_is_unchanged():
    assert remove_quotes_expanded_arg_chars(literal("cat")) == "cat"


def test_double_quoted_name_loses_its_quotes():
    # "cat" must look up the same annotation as cat
    assert remove_quotes_expanded_arg_chars(quoted(literal("cat"))) == "cat"


def test_partially_quoted_name():
    # "c"at -- a string-level strip would produce c"at
    arg_chars = quoted(literal("c")) + literal("at")
    assert remove_quotes_expanded_arg_chars(arg_chars) == "cat"


def test_nested_quotes():
    assert remove_quotes_expanded_arg_chars(quoted(quoted(literal("cat")))) == "cat"


def test_escaped_chars_are_kept():
    arg_chars = [EArgChar(ord("c"))] + literal("at")
    assert remove_quotes_expanded_arg_chars(arg_chars) == "cat"


def test_empty_word():
    assert remove_quotes_expanded_arg_chars([]) == ""
    assert remove_quotes_expanded_arg_chars(quoted([])) == ""


@pytest.mark.parametrize(
    "unexpanded",
    [
        VArgChar("Normal", False, "TOOL", []),
        BArgChar(None),
        AArgChar([]),
    ],
    ids=["variable", "command_substitution", "arithmetic"],
)
def test_unexpanded_words_are_rejected(unexpanded):
    # The value of these is not known at this point, so no name can be
    # produced and the caller must fall back rather than guess.
    assert remove_quotes_expanded_arg_chars([unexpanded]) is None


def test_unexpanded_inside_quotes_is_rejected():
    arg_chars = quoted([VArgChar("Normal", False, "TOOL", [])])
    assert remove_quotes_expanded_arg_chars(arg_chars) is None


def test_unexpanded_alongside_literals_is_rejected():
    arg_chars = literal("pre") + [VArgChar("Normal", False, "X", [])]
    assert remove_quotes_expanded_arg_chars(arg_chars) is None


def test_differs_from_format_arg_chars_only_on_quotes():
    # format_arg_chars regenerates shell syntax and keeps the quotes; this is
    # the behaviour that made the annotation lookup miss for quoted commands.
    arg_chars = quoted(literal("cat"))
    assert format_arg_chars(arg_chars) == '"cat"'
    assert remove_quotes_expanded_arg_chars(arg_chars) == "cat"