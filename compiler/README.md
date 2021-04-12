# The PaSh Compiler

A diagram of the compiler is shown below:

A correspondence between blocks in the diagram and Python modules is shown below:

- PaSh Preprocessor -- [pash.py](../compiler/pash.py)
- Expand, Compile -- [ast_to_ir.py](../compiler/ast_to_ir.py)
- Annotations -- [annotations.py](../compiler/annotations.py), [command_categories.py](../compiler/command_categories.py)
- Optimize -- [pash_runtime.py](../compiler/pash_runtime.py)

**Note:** At the time of the paper submission, PaSh did not have a preprocessing component, and didn't handle variable expansion. These changes significantly improve the practical applicability of PaSh since it can be used on scripts where the environment variables are modified throughout the script.

First, there is the parser in [compiler/parser](../compiler/parser), which is a port of [libdash](https://github.com/mgree/), the dash parser extended with OCaml bindings, extended with ocaml2json and json2ocaml code to interface with PaSh.

Now let's get to the compiler. It's entry point is [compiler/pash.py](../compiler/pash.py) that parses a script and replaces potentially parallelizable regions with calls to [compiler/pash_runtime.sh](../compiler/pash_runtime.sh). It then executes the script.
This allows invoking the compiler during the runtime to have information about the values of environment variables.

The runtime script [compiler/pash_runtime.sh](../compiler/pash_runtime.sh) simply invokes the compiler [compiler/pash_runtime.py](../compiler/pash_runtime.py) and if it succeeds it executes the optimized script, otherwise it executes the original script.

Now the compiler has several stages:

1. It expands words in the AST and then it turns it into our dataflow model (guided by annotations)
   - The expansion and translation happens in [ast_to_ir.py](../compiler/ast_to_ir.py)
   - The dataflow model is mostly defined in [ir.py](../compiler/ir.py)
   - The annotations are processed in [annotations.py](../compiler/annotations.py) and [command_categories.py](../compiler/command_categories.py)
2. It performs transformations on the dataflow graph to expose parallelism (guided by annotations)
   - Translations happen in [pash_runtime.py](../compiler/pash_runtime.py)
3. It then translates the dataflow graph back to a shell script to execute it with bash
   - The `dfg2shell` translation happens in [ir_to_ast.py](../compiler/ir_to_ast.py)
   
 A few interesting fragments are shown below.
 
 The [ast_to_ir.py](https://github.com/andromeda/pash/blob/main/compiler/ast_to_ir.py) contains a case statement that essentially pattern-matches on constructs of the shells script AST and then compiles them accordingly.
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


The following function from [ir.py](https://github.com/andromeda/pash/blob/main/compiler/ir.py) is responsible for parallelizing a single node (i.e., command) in the dataflow graph. Look at the schematic in the comments starting [on line 637](https://github.com/andromeda/pash/blob/main/compiler/ir.py#L637) that gives the high-level overview of what this function does (not shown below).

```Python
    # See comment on line 637
    def parallelize_node(self, node_id, fileIdGen):
        node = self.get_node(node_id)
        assert(node.is_parallelizable())

        ## Initialize the new_node list
        new_nodes = []
        # ... more code ...
```

Another interesting fragment is in [ir_to_ast.py](https://github.com/andromeda/pash/blob/main/compiler/ir_to_ast.py), which translates the parallel dataflow graph back to an AST.

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

