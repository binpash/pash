class UnparallelizableError(Exception):
    pass

class AdjLineNotImplementedError(Exception):
    pass

# to be raised in pash_compiler if a UnparallelizableError is caught at any point running the compiler
#   primarily to differentiate 
#       --assert_compiler_success (exit with error only under general exceptions caught) 
#       --assert_all_regions_parallelizable (exit with error when regions are found not parallelizable + general exceptions)
class NotAllRegionParallelizableError(Exception): 
    pass 