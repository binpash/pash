import sys
import pickle
from ir import *
from json_ast import *

## This file receives the name of a file that holds an IR, reads the
## IR, read some configuration file with node information, and then
## should make a distribution plan for it.

def main():
    # print command line arguments
    ir_filename = sys.argv[1]
    with open(ir_filename, "rb") as ir_file:
        ir_node = pickle.load(ir_file)

    print("Retrieving IR: {} ...".format(ir_filename))
    shell_string = ast_to_shell(ir_node.ast)
    print(shell_string)

    
    ## TODO: As an integration experiment make a distribution planner
    ## that just uses Greenberg's parser to output shell code from the
    ## original AST of the IR, and just execute that.

    ## TODO: In order to be able to execute it, we either have to
    ## execute it in the starting shell (so that we have its state),
    ## or we should somehow pass the parent shell's state to the the
    ## distribution planner, and then the implementation environment.
    ## In general, we probably have to find a way to pass around a
    ## shell's state, as this will be essential for the distributed
    ## setting too.
    
    ## Notes:
    ##
    ## A distribution planner that uses equations (like the ones we
    ## have described with ++) could be separated into three parts:
    ##
    ## 1. Use the equations exhaustively (or until some bound) to
    ##    expose any possible parallelism and distribution. This
    ##    requires that applying equations has a normal form, either
    ##    because they can only be applied a finite number of times,
    ##    or because we have a normal form for equations that can be
    ##    applied an infinite amount of times, that is succinct.
    ##
    ## 2. Assign each node of the IR to a physical node in the system.
    ##
    ## 3. Reapply the equalities (the other way around) so that
    ##    operations on the same node can be merged together to
    ##    improve performance.
        
if __name__ == "__main__":
    main()
    
