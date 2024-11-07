import os
import subprocess
import config

from shell_ast.ast_util import *
from util import ptempfile


##
## This file contains utility functions useful for the speculative execution component
##

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

def good_to_merge(asts):
    return True

def update_df_region(text_to_output: str, trans_option: "SpeculativeTransformationState",
                     df_region_id: int) -> None:
    df_region_path = f'{partial_order_directory()}/{df_region_id}'
    with open(df_region_path, "a") as f:
        f.write(text_to_output)

## TODO: To support partial orders, we need to pass some more context here,
##       i.e., the connections of this node. Now it assumes we have a sequence.
def save_df_region(text_to_output: str, trans_options, df_region_id: int, predecessor_ids: int) -> None:
    ## To support loops we also need to associate nodes with their surrounding loops
    current_bb = trans_options.current_bb()
    log("Df region:", df_region_id, "current bb:", current_bb)

    # Add the loop context to the partial_order state
    trans_options.add_node_bb(df_region_id, current_bb)

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

def serialize_number_of_var_assignments(number_of_var_assignments: int) -> str:
    return f'{number_of_var_assignments}\n'

def serialize_loop_context(node_id: int, bb_id) -> str:
    ## Galaxy brain serialization
    # loop_contexts_str = ",".join([str(loop_ctx) for loop_ctx in loop_contexts])
    bb_id_str = str(bb_id)
    return f'{node_id}-loop_ctx-{bb_id_str}\n'

def serialize_var_assignments(node_id: int) -> str:
    return f'{node_id}-var\n'

def save_current_env_to_file(trans_options):
    initial_env_file = ptempfile()
    subprocess.check_output([f"{os.getenv('PASH_TOP')}/compiler/orchestrator_runtime/pash_declare_vars.sh", initial_env_file])
    partial_order_file_path = trans_options.get_partial_order_file()
    with open(partial_order_file_path, "a") as po_file:
        po_file.write(f'{initial_env_file}\n')

## TODO: Eventually we might want to retrieve the number_of_ids from trans_options
def save_number_of_nodes(trans_options):
    number_of_ids = trans_options.get_number_of_ids()
    partial_order_file_path = trans_options.get_partial_order_file()
    with open(partial_order_file_path, "a") as po_file:
        po_file.write(serialize_number_of_nodes(number_of_ids))

def save_loop_contexts(trans_options):
    node_bb_dict = trans_options.get_all_node_bb()
    log("Loop context dict:", node_bb_dict)
    partial_order_file_path = trans_options.get_partial_order_file()
    with open(partial_order_file_path, "a") as po_file:
        for node_id in sorted(node_bb_dict.keys()):
            bb_id = node_bb_dict[node_id]
            po_file.write(serialize_loop_context(node_id, bb_id))

def save_var_assignment_contexts(trans_options):
    var_nodes = trans_options.get_var_nodes()
    partial_order_file_path = trans_options.get_partial_order_file()
    with open(partial_order_file_path, "a") as po_file:
        for node_id in var_nodes:
            po_file.write(serialize_var_assignments(node_id))
            
def save_number_of_var_assignments(trans_options):
    number_of_var_assignments = trans_options.get_number_of_var_assignments()
    partial_order_file_path = trans_options.get_partial_order_file()
    with open(partial_order_file_path, "a") as po_file:
        po_file.write(serialize_number_of_var_assignments(number_of_var_assignments))

def serialize_partial_order(trans_options):
    ## Initialize the po file
    dir_path = partial_order_directory()
    initialize_po_file(trans_options, dir_path)
    ## Save initial env to po file
    save_current_env_to_file(trans_options)

    ## Save the number of nodes
    save_number_of_nodes(trans_options)

    partial_order_file_path = trans_options.get_partial_order_file()
    with open(partial_order_file_path, "a") as po_file:
        po_file.write('Basic blocks:\n')
        po_file.write('Basic block edges:\n')
        for from_bb_id, to_bb_ids in trans_options.prog.edges.items():
            for to_bb_id, (edge_type, aux_info) in to_bb_ids.items():
                po_file.write(f'{from_bb_id} -> {to_bb_id}: {str(edge_type)}:{aux_info}\n')
        po_file.write('Loop context:\n')

    ## Save loop contexts
    save_loop_contexts(trans_options)
    
    save_number_of_var_assignments(trans_options)
    
    save_var_assignment_contexts(trans_options)

    # Save the edges in the partial order file
    edges = trans_options.get_all_edges()
    with open(partial_order_file_path, "a") as po_file:
        for from_id, to_id in edges:
            po_file.write(serialize_edge(from_id, to_id))
