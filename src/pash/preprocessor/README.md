# PaSh Preprocessor

The PaSh preprocessor transforms shell scripts by identifying candidate dataflow regions and replacing them with calls to the PaSh runtime for compilation and optimization.

## Directory Structure

```
preprocessor/
├── pash_preprocessor.py      # Entry point with all preprocessing logic
├── parse.py                  # Shell script parsing/unparsing
├── util.py                   # Utility functions
├── env_var_names.py          # Environment variable name constants
├── run_tests.sh              # Test runner script
├── shell_ast/
│   ├── ast_to_ast.py              # AST region replacement + preprocess_node
│   ├── walk_preprocess.py         # Generic AST walker with handlers
│   ├── transformation_options.py  # Transformation state classes
│   ├── ast_util.py                # AST utility functions
│   ├── handlers/
│   │   └── loop_tracking.py       # ForNode loop tracking handler
│   └── test_walk_preprocess.py    # Unit tests
└── speculative/
    └── util_spec.py               # Speculative execution support
```

## Running the Preprocessor

### Via pash command (normal usage)

```bash
pash script.sh
```

### Standalone preprocessing

```bash
export PASH_TOP=/path/to/pash
export PASH_TMP_PREFIX=$(mktemp -d)/
export PASH_BASH_VERSION="5 2 21"
export PYTHONPATH="$PASH_TOP/preprocessor:$PASH_TOP/compiler:$PYTHONPATH"

python3 preprocessor/pash_preprocessor.py --output /tmp/out.sh script.sh
```

## Running Tests

### Using the test script

```bash
cd preprocessor
./run_tests.sh
```

### With verbose output

```bash
./run_tests.sh -v
```

### Run specific test class

```bash
./run_tests.sh TestPreprocessNode
```

### Run specific test method

```bash
./run_tests.sh TestPreprocessNode.test_pipe_node
```

### Manual test execution

```bash
export PASH_TOP=/path/to/pash
export PASH_TMP_PREFIX=$(mktemp -d)/
export PASH_BASH_VERSION="5 2 21"
export PYTHONPATH="$PASH_TOP/preprocessor:$PASH_TOP/compiler:$PYTHONPATH"

cd preprocessor
python3 -m unittest shell_ast.test_walk_preprocess -v
```

## Architecture

### Preprocessing Flow

1. **Entry Point** (`pash_preprocessor.py`): Parses arguments, orchestrates parsing, preprocessing, and unparsing
2. **Parsing** (`parse.py`): Converts shell script to AST using libdash/libbash
3. **AST Walking** (`walk_preprocess.py`): Generic pattern-matching walker with handler support
4. **Region Replacement** (`ast_to_ast.py`): Identifies dataflow regions and replaces with runtime calls
5. **Unparsing** (`parse.py`): Converts transformed AST back to shell script

### Key Classes

- `WalkPreprocess`: Generic AST walker using Python pattern matching
- `PreprocessContext`: Context threaded through traversal (trans_options, last_object)
- `NodeResult`: Result from processing a node (ast, replace_whole, non_maximal, something_replaced)
- `PreprocessedAST`: Result wrapper with flags for replacement decisions
- `TransformationState`: Manages node IDs, loop contexts, and replacement generation

### Node Categories

| Category | Nodes | Behavior |
|----------|-------|----------|
| Leaf Replacements | `PipeNode`, `CommandNode`, `BackgroundNode` | Marked for replacement |
| Binary Operators | `SemiNode`, `AndNode`, `OrNode` | Close-recurse both children |
| Single Child | `RedirNode`, `SubshellNode`, `NotNode`, `GroupNode`, etc. | Close-recurse single child |
| Control Flow | `IfNode`, `WhileNode`, `ForNode`, `CaseNode` | Close-recurse with special handling |
| No-Op | `DefunNode`, `ArithNode` | Return unchanged |

### Custom Handlers

The walker supports custom handlers for specific node types. For example, `ForNode` has a custom handler in `handlers/loop_tracking.py` that injects loop iteration tracking code for speculative execution support.

```python
from shell_ast.walk_preprocess import WalkPreprocess

handlers = {"for": my_custom_for_handler}
walker = WalkPreprocess(handlers=handlers)
```

## Dependencies

- `shasta`: AST node definitions
- `libdash`: POSIX shell parser
- `libbash`: Bash parser (for `--bash` mode)

## Environment Variables

- `PASH_TOP`: Root directory of PaSh installation
- `PASH_TMP_PREFIX`: Temporary directory for preprocessor output
- `PASH_BASH_VERSION`: Bash version string (e.g., "5 2 21")
