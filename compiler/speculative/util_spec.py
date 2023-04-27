
import os
import config

from shell_ast.ast_util import *

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
    
    def get_current_id(self):
        return self.counter - 1
    
    def get_number_of_ids(self):
        return self.counter

## TODO: Move this to the trans_options class
ID_GENERATOR = IdGen()

def initialize(trans_options) -> None:
    ## Make the directory that contains the files in the partial order
    dir_path = partial_order_directory()
    os.makedirs(dir_path)
    # ## Initialize the po file
    # initialize_po_file(trans_options, dir_path)

def partial_order_directory() -> str:
    return f'{config.PASH_TMP_PREFIX}/speculative/partial_order/'

def partial_order_file_path():
    return f'{config.PASH_TMP_PREFIX}/speculative/partial_order_file'

def initialize_po_file(trans_options, dir_path) -> None:
    ## Initializae the partial order file
    with open(trans_options.get_partial_order_file(), 'w') as f:
        f.write(f'# Partial order files path:\n')
        f.write(f'{dir_path}\n')

def scheduler_server_init_po_msg(partial_order_file: str) -> str:
    return f'Init:{partial_order_file}'

def get_next_id():
    global ID_GENERATOR
    return ID_GENERATOR.get_next_id()

## TODO: To support partial orders, we need to pass some more context here,
##       i.e., the connections of this node. Now it assumes we have a sequence.
def save_df_region(text_to_output: str, trans_options, df_region_id: int, predecessor_ids: int) -> None:
    ## To support loops we also need to associate nodes with their surrounding loops
    current_loop_context = trans_options.get_current_loop_context()
    log("Df region:", df_region_id, "loop context:", current_loop_context)

    # Add the loop context to the partial_order state
    trans_options.add_node_loop_context(df_region_id, current_loop_context)

    # Save df_region as text in its own file
    df_region_path = f'{partial_order_directory()}/{df_region_id}'
    with open(df_region_path, "w") as f:
        f.write(text_to_output)

    ## Save the edges in the partial order state
    for predecessor in predecessor_ids:
        trans_options.add_edge(predecessor, df_region_id)



## TODO: Figure out a way to put all serialization/deserialization of messages
##       and parsing/unparsing in a specific module.

## TODO: Move serialization to a partial_order_file.py
def serialize_edge(from_id: int, to_id: int) -> str:
    return f'{from_id} -> {to_id}\n'

def serialize_number_of_nodes(number_of_ids: int) -> str:
    return f'{number_of_ids}\n'

## TODO: Eventually we might want to retrieve the number_of_ids from trans_options
def save_number_of_nodes(trans_options):
    number_of_ids = ID_GENERATOR.get_number_of_ids()
    partial_order_file_path = trans_options.get_partial_order_file()
    with open(partial_order_file_path, "a") as po_file:
        po_file.write(serialize_number_of_nodes(number_of_ids))

def serialize_partial_order(trans_options):
    ## Initialize the po file
    dir_path = partial_order_directory()
    initialize_po_file(trans_options, dir_path)

    # Save the edges in the partial order file
    partial_order_file_path = trans_options.get_partial_order_file()
    edges = trans_options.get_all_edges()
    with open(partial_order_file_path, "a") as po_file:
        for from_id, to_id in edges:
            po_file.write(serialize_edge(from_id, to_id))

    ## Save the number of nodes
    ## TODO: Move that in the beginning of the file
    save_number_of_nodes(trans_options)
