
import os
import config

##
## This file contains utility functions useful for the speculative execution component
##

## TODO: There is a similar class in ir.py. Could we avoid duplication?
class IdGen:
    def __init__(self, counter=0):
        self.counter = counter
    
    def get_next_id(self):
        new_id = self.counter
        self.counter += 1
        return new_id

## TODO: Should we move this to the trans_options class 
##       (which we could rename to trans_config) and make a subclass for
##       the two different transformations.
ID_GENERATOR = IdGen()

def initialize(trans_options) -> None:
    ## Make the directory that contains the files in the partial order
    dir_path = partial_order_directory()
    os.makedirs(dir_path)
    ## Initialize the po file
    initialize_po_file(trans_options, dir_path)

def partial_order_directory() -> str:
    return f'{config.PASH_TMP_PREFIX}/speculative/partial_order/'

def partial_order_file_path():
    return f'{config.PASH_TMP_PREFIX}/speculative/partial_order_file'

def initialize_po_file(trans_options, dir_path) -> None:
    ## Initializae the partial order file
    with open(trans_options.get_partial_order_file(), 'w') as f:
        f.write(f'# Partial order files path:\n')
        f.write(f'{dir_path}\n')

def get_next_id():
    global ID_GENERATOR
    return ID_GENERATOR.get_next_id()

## TODO: To support partial orders, we need to pass some more context here,
##       i.e., the connections of this node. Now it assumes we have a sequence.
def save_df_region(text_to_output: str, trans_options, df_region_id: int, predecessor_ids: int) -> None:
    # Save df_region as text in its own file
    df_region_path = f'{partial_order_directory()}/{df_region_id}'
    with open(df_region_path, "w") as f:
        f.write(text_to_output)

    # Save the edges in the partial order file
    partial_order_file_path = trans_options.get_partial_order_file()
    with open(partial_order_file_path, "a") as po_file:
        for predecessor in predecessor_ids:
            po_file.write(serialize_edge(predecessor, df_region_id))

def serialize_edge(from_id, to_id):
    return f'{from_id} -> {to_id}\n'
