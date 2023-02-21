# The PaSh Compiler
Quick Jump: [intro](#introduction) | [overview](#compiler-overview) | [details](#zooming-into-fragments) | [earlier versions](#earlier-versions)

## Introduction

PaSh has recently shifted away from ahead-of-time compilation and towards just-in-time compilation intermixed with the execution of a script.
This shift brings many benefits, allowing PaSh to correctly handle expansion and other important details -- but complicates the clear exposition of the two phases.
A high-level diagram of PaSh's end-to-end operation is shown below:

<img src="https://docs.google.com/drawings/d/e/2PACX-1vSIuacgBR_QFOzawoAdJMmTjgsdnDUkp1DbSjLVlrowlhL6kxqckXXsL7SPoRXKfaC1hw9HQzJitmDP/pub?w=1364&amp;h=454">

PaSh pre-processes a sequential script to insert calls to the `pash_compiler.py`.
It then invokes the script, switching between evaluation, execution, and parallelization at runtime:
(i) it first parses the script, creating an abstact syntax tree (AST); 
(ii) it then expands the nodes of the AST, often calling the shell which performs that expansion;
(iii) it compiles dataflow regions, parts of the AST that are potentially parallelizable, through an iterative optimization proceedure applied over a dataflow graph (DFG); and
(iv) finally emits the parallel script by translating the DFG to AST and unparsing the AST back to a shell script.
The compilation takes into account information about individual commands through annotations, and the emitted parallel script uses additional constructs provided by PaSh's [runtime library](../runtime).

A correspondence between blocks in the diagram and Python modules is shown below:

- Preprocessing: [pash.py](./pash.py)
- Expansion and compilation: [ast_to_ir.py](./ast_to_ir.py)
- Optimization: [pash_compiler.py](./pash_compiler.py)

## Compiler Overview

Now let's get to the compiler.
It's entry point is [pash.py](./pash.py) that parses a script and replaces potentially parallelizable regions with calls to [pash_runtime.sh](./pash_runtime.sh).
It then executes the script.
This allows invoking the compiler during the runtime to have information about the values of environment variables.

The [pash_runtime.sh](./pash_runtime.sh) script simply invokes the [pash.py](./pash.py) compiler:
  if it succeeds it executes the optimized script, otherwise it executes the original script.

The compiler has several stages:

1. It expands words in the AST and then it turns it into our dataflow model (guided by annotations)
   - The expansion and translation happens in [ast_to_ir.py](./ast_to_ir.py)
   - The dataflow model is defined mostly in [ir.py](./ir.py)
   - The annotations are processed in [binpash/annotations](https://github.com/binpash/annotations)
2. It performs transformations on the dataflow graph to expose parallelism (guided by annotations)
   - Translations happen in [pash_compiler.py](./pash_compiler.py)
3. It then translates the dataflow graph back to a shell script to execute it with bash
   - The `dfg2shell` translation happens in [ir_to_ast.py](./ir_to_ast.py)

[//]: # (TODO: the parsing/unparsing components need update)

## Zooming into Fragments
   
A few interesting fragments are outlined below.
 
The [ast_to_ir.py](./ast_to_ir.py) file contains a case statement that essentially pattern-matches on constructs of the shells script AST and then compiles them accordingly.
```Python
 compile_cases = {
        "Pipe": (lambda fileIdGen, config:
                 lambda ast_node: compile_node_pipe(ast_node, fileIdGen, config)),
        "Command": (lambda fileIdGen, config:
                    lambda ast_node: compile_node_command(ast_node, fileIdGen, config)),
        "And": (lambda fileIdGen, config:
                lambda ast_node: compile_node_and_or_semi(ast_node, fileIdGen, config)),
        "Or": (lambda fileIdGen, config:
               lambda ast_node: compile_node_and_or_semi(ast_node, fileIdGen, config))
        # ... more code ...
```

The following function from [ir.py](./ir.py) is responsible for parallelizing a single node (_i.e._, a command) in the dataflow graph.
Look at the schematic in the comments starting [on line 663](./ir.py#L663) that gives the high-level overview of what this function does (not shown below).

[//]: # (TODO: Add schematic here)

```Python
    # See comment on line 637
    def parallelize_node(self, node_id, fileIdGen):
        node = self.get_node(node_id)
        assert(node.is_parallelizable())

        ## Initialize the new_node list
        new_nodes = []
        # ... more code ...
```

Another interesting fragment is in [ir_to_ast.py](./ir_to_ast.py), which translates the parallel dataflow graph back to an AST.

```Python
def ir2ast(ir, args):
    # ... more code ...
    
    ## Make the main body
    body = ir.to_ast(drain_streams)
    
    # ... more code ...
    
    ## Call the prologue that creates fifos for all ephemeral fids    
    prologue = make_ir_prologue(ephemeral_fids)
    
    ## Call the epilogue that removes all ephemeral fids
    epilogue = make_ir_epilogue(ephemeral_fids, clean_up_graph, args.log_file)

    final_asts = prologue + body + epilogue
```

This AST is then unparsed back into a (parallel) shell script.

## Earlier Versions

The compiler is outlined in the [EuroSys paper](https://arxiv.org/pdf/2007.09436.pdf), but has evolved considerably since then:

* PaSh originally did not have a preprocessing component, and didn't handle variable expansion. It now does both, significantly improving its practical applicability since it can be used on scripts where the environment variables are modified throughout the script.

* PaSh originally was using code in [parser](./parser) -- a port of [libdash](https://github.com/mgree/), the `dash` parser extended with OCaml bindings -- and specifically the `ocaml2json` and `json2ocaml` binaries to interface with PaSh. PaSh now uses a custom parser written in Python, avoiding any dependency to OCaml and simplifying dependency management.
 
