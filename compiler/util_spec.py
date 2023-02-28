
##
## This file contains utility functions useful for the speculative execution component
##

def initialize_po_file(trans_options) -> None:
    ## Initializae the partial order file
    open(trans_options.get_partial_order_file(), 'w').close()

