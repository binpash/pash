"""
Test cases for walk_preprocess module.

This module tests the refactored preprocessing walker against various
shell script constructs to ensure correct behavior.
"""

import os
import sys
import tempfile
import unittest

# Add compiler directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Set required environment variables before importing modules
os.environ.setdefault("PASH_TOP", os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))))
os.environ.setdefault("PASH_TMP_PREFIX", tempfile.mkdtemp() + "/")
os.environ.setdefault("PASH_BASH_VERSION", "5 2 21")

from parse import parse_shell_to_asts
from shell_ast.transformation_options import TransformationState
from shell_ast.preprocess_ast_cases import preprocess_node, preprocess_close_node
from preprocessor.preprocessor import preprocess


class MockArgs:
    """Mock arguments for preprocessing."""
    def __init__(self):
        self.bash = False
        self.preprocess_mode = "pash"


def parse_script(script: str):
    """Helper to parse a shell script string to AST objects."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.sh', delete=False) as f:
        f.write(script)
        f.flush()
        script_path = f.name

    try:
        return list(parse_shell_to_asts(script_path, bash_mode=False))
    finally:
        os.unlink(script_path)


def preprocess_script(script: str) -> str:
    """Helper to preprocess a script and return the result."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.sh', delete=False) as f:
        f.write(script)
        f.flush()
        script_path = f.name

    try:
        args = MockArgs()
        return preprocess(script_path, args)
    finally:
        os.unlink(script_path)


class TestPreprocessNode(unittest.TestCase):
    """Test preprocess_node function for individual node types."""

    def test_command_node(self):
        """Test that CommandNode with arguments is marked for replacement."""
        ast_objects = parse_script("echo hello")
        self.assertEqual(len(ast_objects), 1)

        ast, _, _, _ = ast_objects[0]
        trans_options = TransformationState()
        result = preprocess_node(ast, trans_options, last_object=True)

        self.assertEqual(result.replace_whole, True)
        self.assertEqual(result.non_maximal, False)
        self.assertEqual(result.something_replaced, True)

    def test_command_assignment_only(self):
        """Test that assignment-only CommandNode is not marked for replacement."""
        ast_objects = parse_script("x=1")
        self.assertEqual(len(ast_objects), 1)

        ast, _, _, _ = ast_objects[0]
        trans_options = TransformationState()
        result = preprocess_node(ast, trans_options, last_object=True)

        self.assertEqual(result.replace_whole, False)
        self.assertEqual(result.something_replaced, False)

    def test_pipe_node(self):
        """Test that PipeNode is marked for replacement."""
        ast_objects = parse_script("echo hello | grep h")
        self.assertEqual(len(ast_objects), 1)

        ast, _, _, _ = ast_objects[0]
        trans_options = TransformationState()
        result = preprocess_node(ast, trans_options, last_object=True)

        self.assertEqual(result.replace_whole, True)
        self.assertEqual(result.non_maximal, False)
        self.assertEqual(result.something_replaced, True)

    def test_background_node(self):
        """Test that BackgroundNode is marked as non-maximal."""
        ast_objects = parse_script("cmd &")
        self.assertEqual(len(ast_objects), 1)

        ast, _, _, _ = ast_objects[0]
        trans_options = TransformationState()
        result = preprocess_node(ast, trans_options, last_object=True)

        self.assertEqual(result.replace_whole, True)
        self.assertEqual(result.non_maximal, True)
        self.assertEqual(result.something_replaced, True)

    def test_semi_node(self):
        """Test that SemiNode recurses into children."""
        ast_objects = parse_script("echo hello; echo world")
        self.assertEqual(len(ast_objects), 1)

        ast, _, _, _ = ast_objects[0]
        trans_options = TransformationState()
        result = preprocess_node(ast, trans_options, last_object=True)

        self.assertEqual(result.replace_whole, False)
        self.assertEqual(result.something_replaced, True)

    def test_and_node(self):
        """Test that AndNode recurses into children."""
        ast_objects = parse_script("cmd1 && cmd2")
        self.assertEqual(len(ast_objects), 1)

        ast, _, _, _ = ast_objects[0]
        trans_options = TransformationState()
        result = preprocess_node(ast, trans_options, last_object=True)

        self.assertEqual(result.replace_whole, False)
        self.assertEqual(result.something_replaced, True)

    def test_or_node(self):
        """Test that OrNode recurses into children."""
        ast_objects = parse_script("cmd1 || cmd2")
        self.assertEqual(len(ast_objects), 1)

        ast, _, _, _ = ast_objects[0]
        trans_options = TransformationState()
        result = preprocess_node(ast, trans_options, last_object=True)

        self.assertEqual(result.replace_whole, False)
        self.assertEqual(result.something_replaced, True)

    def test_not_node(self):
        """Test that NotNode recurses into body."""
        ast_objects = parse_script("! cmd")
        self.assertEqual(len(ast_objects), 1)

        ast, _, _, _ = ast_objects[0]
        trans_options = TransformationState()
        result = preprocess_node(ast, trans_options, last_object=True)

        self.assertEqual(result.replace_whole, False)
        self.assertEqual(result.something_replaced, True)

    def test_if_node(self):
        """Test that IfNode recurses into all branches."""
        ast_objects = parse_script("if true; then echo yes; else echo no; fi")
        self.assertEqual(len(ast_objects), 1)

        ast, _, _, _ = ast_objects[0]
        trans_options = TransformationState()
        result = preprocess_node(ast, trans_options, last_object=True)

        self.assertEqual(result.replace_whole, False)
        self.assertEqual(result.something_replaced, True)

    def test_for_node(self):
        """Test that ForNode recurses into body."""
        ast_objects = parse_script("for x in 1 2 3; do echo $x; done")
        self.assertEqual(len(ast_objects), 1)

        ast, _, _, _ = ast_objects[0]
        trans_options = TransformationState()
        result = preprocess_node(ast, trans_options, last_object=True)

        self.assertEqual(result.replace_whole, False)
        self.assertEqual(result.something_replaced, True)

    def test_while_node(self):
        """Test that WhileNode recurses into test and body."""
        ast_objects = parse_script("while true; do echo x; done")
        self.assertEqual(len(ast_objects), 1)

        ast, _, _, _ = ast_objects[0]
        trans_options = TransformationState()
        result = preprocess_node(ast, trans_options, last_object=True)

        self.assertEqual(result.replace_whole, False)
        self.assertEqual(result.something_replaced, True)

    def test_subshell_node(self):
        """Test that SubshellNode recurses into body."""
        ast_objects = parse_script("(echo subshell)")
        self.assertEqual(len(ast_objects), 1)

        ast, _, _, _ = ast_objects[0]
        trans_options = TransformationState()
        result = preprocess_node(ast, trans_options, last_object=True)

        self.assertEqual(result.replace_whole, False)
        self.assertEqual(result.something_replaced, True)

    def test_group_node_as_command(self):
        """Test brace group parsing behavior.

        Note: The POSIX (dash) parser optimizes away group nodes for simple
        cases like '{ echo group; }', turning them into plain CommandNode.
        GroupNode is only produced by the bash parser.
        This test verifies we handle what we get from the parser.
        """
        # With dash parser, '{ echo group; }' parses as CommandNode
        ast_objects = parse_script("{ echo group; }")
        self.assertEqual(len(ast_objects), 1)

        ast, _, _, _ = ast_objects[0]
        trans_options = TransformationState()
        result = preprocess_node(ast, trans_options, last_object=True)

        # The parser returns CommandNode, which should be replaced
        self.assertEqual(type(ast).__name__, "CommandNode")
        self.assertEqual(result.replace_whole, True)
        self.assertEqual(result.something_replaced, True)

    def test_case_node(self):
        """Test that CaseNode recurses into case bodies."""
        ast_objects = parse_script("case $x in a) echo a;; esac")
        self.assertEqual(len(ast_objects), 1)

        ast, _, _, _ = ast_objects[0]
        trans_options = TransformationState()
        result = preprocess_node(ast, trans_options, last_object=True)

        self.assertEqual(result.replace_whole, False)
        self.assertEqual(result.something_replaced, True)

    def test_defun_node(self):
        """Test that DefunNode is not preprocessed (per current logic)."""
        ast_objects = parse_script("myfunc() { echo hello; }")
        self.assertEqual(len(ast_objects), 1)

        ast, _, _, _ = ast_objects[0]
        trans_options = TransformationState()
        result = preprocess_node(ast, trans_options, last_object=True)

        self.assertEqual(result.replace_whole, False)
        self.assertEqual(result.something_replaced, False)


class TestFullPreprocessing(unittest.TestCase):
    """Test full preprocessing pipeline."""

    def test_simple_command(self):
        """Test preprocessing of simple command."""
        result = preprocess_script("echo hello")
        self.assertIn("pash_runtime", result)

    def test_pipeline(self):
        """Test preprocessing of pipeline."""
        result = preprocess_script("echo hello | grep h")
        self.assertIn("pash_runtime", result)

    def test_for_loop_injects_tracking(self):
        """Test that for loop preprocessing injects loop tracking."""
        result = preprocess_script("for x in 1 2 3; do echo $x; done")
        self.assertIn("pash_loop_0_iter", result)
        self.assertIn("pash_loop_iters", result)

    def test_nested_for_loops(self):
        """Test that nested for loops have different loop IDs."""
        result = preprocess_script("for x in 1 2; do for y in a b; do echo $x $y; done; done")
        self.assertIn("pash_loop_0_iter", result)
        self.assertIn("pash_loop_1_iter", result)

    def test_assignment_preserved(self):
        """Test that assignment without command is preserved as-is."""
        result = preprocess_script("x=1")
        self.assertEqual(result.strip(), "x=1")

    def test_if_statement(self):
        """Test preprocessing of if statement."""
        result = preprocess_script("if true; then echo yes; fi")
        self.assertIn("pash_runtime", result)

    def test_complex_script(self):
        """Test preprocessing of complex script."""
        script = """
echo start
for x in 1 2 3; do
    echo $x | grep 1
done
if true; then
    cat file | sort
else
    cat file | uniq
fi
echo end
"""
        result = preprocess_script(script)
        self.assertIn("pash_runtime", result)
        self.assertIn("pash_loop_0_iter", result)


class TestPreprocessCloseNode(unittest.TestCase):
    """Test preprocess_close_node function."""

    def test_close_node_replaces_command(self):
        """Test that close_node replaces command nodes.

        Note: The replacement is still a CommandNode, but it's a different
        CommandNode that sources the pash_runtime.sh script.
        """
        ast_objects = parse_script("echo hello")
        ast, _, _, _ = ast_objects[0]
        trans_options = TransformationState()

        final_ast, something_replaced = preprocess_close_node(
            ast, trans_options, last_object=True
        )

        self.assertTrue(something_replaced)
        # The returned AST is still a CommandNode but with different content
        # (it's the runtime call). We verify this by checking it's different
        # from the original.
        self.assertIsNot(final_ast, ast)

    def test_close_node_does_not_replace_assignment(self):
        """Test that close_node does not replace assignment-only nodes."""
        ast_objects = parse_script("x=1")
        ast, _, _, _ = ast_objects[0]
        trans_options = TransformationState()

        final_ast, something_replaced = preprocess_close_node(
            ast, trans_options, last_object=True
        )

        self.assertFalse(something_replaced)
        # The AST should be unchanged
        self.assertEqual(type(final_ast).__name__, "CommandNode")


if __name__ == "__main__":
    unittest.main()
