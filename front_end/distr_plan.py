import sys
import pickle
from ir import *

## This file receives the name of a file that holds an IR, reads the
## IR, read some configuration file with node information, and then
## should make a distribution plan for it.

def main():
    # print command line arguments
    ir_filename = sys.argv[1]
    with open(ir_filename, "rb") as ir_file:
        ir_node = pickle.load(ir_file)

    print(ir_node)

    ## TODO: Do something with the node
    
        
if __name__ == "__main__":
    main()
    
