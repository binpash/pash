# PaSh Preprocessor

The PaSh preprocessor transforms shell scripts by identifying candidate dataflow regions and replacing them with calls to the PaSh runtime for compilation and optimization.

## Directory Structure

```
preprocessor/
├── pash_preprocessor.py          # Entry point + region replacement logic
├── parse.py                      # Shell script parsing/unparsing
├── util.py                       # Utility functions
├── ast_util.py                   # AST construction helpers
├── transformation_options.py     # Transformation state classes
├── walk_preprocess.py            # Preprocessing walker (extends shasta's CommandVisitor)
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

## Architecture

### Preprocessing Flow

1. **Entry Point** (`pash_preprocessor.py`): Parses arguments, orchestrates parsing, preprocessing, and unparsing
2. **Parsing** (`parse.py`): Converts shell script to AST using libdash/libbash
3. **AST Walking** (`walk_preprocess.py`): Preprocessing visitor built on `shasta.ast_walker.CommandVisitor`
4. **Region Replacement** (`pash_preprocessor.py`): Identifies dataflow regions and replaces with runtime calls
5. **Unparsing** (`parse.py`): Converts transformed AST back to shell script

### Key Classes

- `CommandVisitor` (from `shasta`): Generic command-level AST visitor base class
- `WalkPreprocess`: Preprocessing visitor extending `CommandVisitor` with close-node semantics
- `PreprocessContext`: Context threaded through traversal (trans_options, last_object)
- `NodeResult`: Result from processing a node (ast, replace_whole, non_maximal, something_replaced)
- `PreprocessedAST`: Result wrapper with flags for replacement decisions
- `TransformationState`: Manages node IDs and replacement generation

### Node Categories

| Category | Nodes | Behavior |
|----------|-------|----------|
| Leaf Replacements | `PipeNode`, `CommandNode`, `BackgroundNode` | Marked for replacement |
| Binary Operators | `SemiNode`, `AndNode`, `OrNode` | Close-recurse both children |
| Single Child | `RedirNode`, `SubshellNode`, `NotNode`, `GroupNode`, etc. | Close-recurse single child |
| Control Flow | `IfNode`, `WhileNode`, `ForNode`, `CaseNode` | Close-recurse all children |
| No-Op | `DefunNode`, `ArithNode` | Return unchanged |

## Dependencies

- `shasta`: AST node definitions and generic `CommandVisitor`
- `libdash`: POSIX shell parser
- `libbash`: Bash parser (for `--bash` mode)

## Environment Variables

- `PASH_TOP`: Root directory of PaSh installation
- `PASH_TMP_PREFIX`: Temporary directory for preprocessor output
- `PASH_BASH_VERSION`: Bash version string (e.g., "5 2 21")
