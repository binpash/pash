import os

from ir import *
from util import *

## Ensure that PASH_TMP_PREFIX is set by pa.sh
assert not os.getenv("PASH_TIMESTAMP") is None
PASH_TIMESTAMP = os.getenv("PASH_TIMESTAMP")
DIR_NAME = f"pash_graphviz_{PASH_TIMESTAMP}"


def maybe_init_graphviz_dir(args):
    if not args.graphviz == "no":
        init_graphviz_dir(args)


def init_graphviz_dir(args):
    graphviz_dir_path = os.path.join(args.graphviz_dir, DIR_NAME)

    try:
        os.mkdir(graphviz_dir_path)
    except:
        print(f"Error: Graphviz dir:{graphviz_dir_path} could not be created!")
        exit(1)

    log("Created graphviz dir:", graphviz_dir_path)


def maybe_generate_graphviz(ir: IR, args, name="dfg"):
    if not args.graphviz == "no":
        generate_graphviz(ir, args, name=name)


def generate_graphviz(ir: IR, args, name="dfg"):
    ## TODO: It is unclear if importing in here (instead of in general)
    ##       improves startup cost of the pash_runtime when not using graphviz.
    import graphviz

    log("Generating graph...")
    dot = ir.generate_graphviz()

    ## The option argument of graphviz contains the format
    dot.format = args.graphviz

    graphviz_dir_path = os.path.join(args.graphviz_dir, DIR_NAME)
    dot.render(directory=graphviz_dir_path, filename=name)

    log("Saved graph visualization in:", os.path.join(graphviz_dir_path, name))
