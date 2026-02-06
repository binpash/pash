"""
Utility functions for the speculative execution component.

Adapted from spec_future branch to work with main branch's module structure,
while producing the partial order file format expected by hs scheduler_server.py.
"""

import os
import subprocess

from util import log, PASH_TMP_PREFIX, ptempfile


def initialize(trans_options) -> None:
    """Initialize the partial order directory."""
    dir_path = partial_order_directory()
    os.makedirs(dir_path)


def partial_order_directory() -> str:
    """Return the path to the partial order directory."""
    return f"{PASH_TMP_PREFIX}/speculative/partial_order/"


def partial_order_file_path():
    """Return the path to the partial order file."""
    return f"{PASH_TMP_PREFIX}/speculative/partial_order_file"


def initialize_po_file(trans_options, dir_path) -> None:
    """Initialize the partial order file."""
    with open(trans_options.get_partial_order_file(), "w") as f:
        f.write(f"# Partial order files path:\n")
        f.write(f"{dir_path}\n")


def scheduler_server_init_po_msg(partial_order_file: str) -> str:
    """Create message to initialize scheduler with partial order file."""
    return f"Init:{partial_order_file}"


def save_df_region(
    text_to_output: str, trans_options, df_region_id: int, predecessor_ids: int
) -> None:
    """Save a dataflow region to a file."""
    # Associate nodes with their surrounding loops
    current_loop_context = trans_options.get_current_loop_context()
    log("Df region:", df_region_id, "loop context:", current_loop_context)

    # Add the loop context to the partial_order state
    trans_options.add_node_loop_context(df_region_id, current_loop_context)

    # Save df_region as text in its own file
    df_region_path = f"{partial_order_directory()}/{df_region_id}"
    with open(df_region_path, "w", encoding="utf-8") as f:
        f.write(text_to_output)

    # Save the edges in the partial order state
    for predecessor in predecessor_ids:
        trans_options.add_edge(predecessor, df_region_id)


def serialize_edge(from_id: int, to_id: int) -> str:
    """Serialize an edge in the partial order."""
    return f"{from_id} -> {to_id}\n"


def serialize_number_of_nodes(number_of_ids: int) -> str:
    """Serialize the number of nodes."""
    return f"{number_of_ids}\n"


def serialize_number_of_var_assignments(number_of_var_assignments: int) -> str:
    """Serialize the number of variable assignments."""
    return f"{number_of_var_assignments}\n"


def serialize_loop_context(node_id: int, loop_contexts) -> str:
    """Serialize a loop context."""
    loop_contexts_str = ",".join([str(loop_ctx) for loop_ctx in loop_contexts])
    return f"{node_id}-loop_ctx-{loop_contexts_str}\n"


def serialize_var_assignments(node_id: int) -> str:
    """Serialize a variable assignment node."""
    return f"{node_id}-var\n"


def save_current_env_to_file(trans_options):
    """Save the current environment to a file and record it in the partial order."""
    initial_env_file = ptempfile()
    ## Use PASH_SPEC_TOP for the declare_vars script (in the hs repo's jit_runtime)
    pash_spec_top = os.getenv('PASH_SPEC_TOP', '')
    declare_vars_script = os.path.join(pash_spec_top, 'jit_runtime', 'pash_declare_vars.sh')
    if os.path.exists(declare_vars_script):
        subprocess.check_output([declare_vars_script, initial_env_file])
    else:
        ## Fallback: create a minimal env file
        log("Warning: pash_declare_vars.sh not found at", declare_vars_script)
        with open(initial_env_file, 'w') as f:
            f.write("")
    partial_order_file_path = trans_options.get_partial_order_file()
    with open(partial_order_file_path, "a") as po_file:
        po_file.write(f"{initial_env_file}\n")


def save_number_of_nodes(trans_options):
    """Save the number of nodes to the partial order file."""
    number_of_ids = trans_options.get_number_of_ids()
    po_file_path = trans_options.get_partial_order_file()
    with open(po_file_path, "a") as po_file:
        po_file.write(serialize_number_of_nodes(number_of_ids))


def save_loop_contexts(trans_options):
    """Save loop contexts to the partial order file."""
    loop_context_dict = trans_options.get_all_loop_contexts()
    log("Loop context dict:", loop_context_dict)
    po_file_path = trans_options.get_partial_order_file()
    with open(po_file_path, "a") as po_file:
        for node_id in sorted(loop_context_dict.keys()):
            loop_ctx = loop_context_dict[node_id]
            po_file.write(serialize_loop_context(node_id, loop_ctx))


def save_var_assignment_contexts(trans_options):
    """Save variable assignment contexts to the partial order file."""
    if hasattr(trans_options, 'get_var_nodes'):
        var_nodes = trans_options.get_var_nodes()
        po_file_path = trans_options.get_partial_order_file()
        with open(po_file_path, "a") as po_file:
            for node_id in var_nodes:
                po_file.write(serialize_var_assignments(node_id))


def save_number_of_var_assignments(trans_options):
    """Save number of variable assignments to the partial order file."""
    if hasattr(trans_options, 'get_number_of_var_assignments'):
        number_of_var_assignments = trans_options.get_number_of_var_assignments()
    else:
        number_of_var_assignments = 0
    po_file_path = trans_options.get_partial_order_file()
    with open(po_file_path, "a") as po_file:
        po_file.write(serialize_number_of_var_assignments(number_of_var_assignments))


def serialize_partial_order(trans_options):
    """Serialize the complete partial order to a file.

    Format expected by hs scheduler_server.py:
    1. cmds_directory
    2. initial_env_file
    3. number_of_nodes
    4. "Basic blocks:" header
    5. "Basic block edges:" header + edges
    6. "Loop context:" header + contexts
    7. number_of_var_assignments
    8. var assignments
    9. edges
    """
    # Initialize the po file (writes directory path)
    dir_path = partial_order_directory()
    initialize_po_file(trans_options, dir_path)

    # Save initial env to po file
    save_current_env_to_file(trans_options)

    # Save the number of nodes
    save_number_of_nodes(trans_options)

    po_file_path = trans_options.get_partial_order_file()
    with open(po_file_path, "a") as po_file:
        po_file.write("Basic blocks:\n")
        po_file.write("Basic block edges:\n")

        # Write basic block edges if available
        if hasattr(trans_options, 'prog') and hasattr(trans_options.prog, 'edges'):
            for from_bb_id, to_bb_ids in trans_options.prog.edges.items():
                for to_bb_id, edge_type in to_bb_ids.items():
                    po_file.write(f"{from_bb_id} -> {to_bb_id}: {str(edge_type)}\n")

        po_file.write("Loop context:\n")

    # Save loop contexts
    save_loop_contexts(trans_options)

    # Save var assignments
    save_number_of_var_assignments(trans_options)
    save_var_assignment_contexts(trans_options)

    # Save the edges in the partial order file
    edges = trans_options.get_all_edges()
    with open(po_file_path, "a") as po_file:
        for from_id, to_id in edges:
            po_file.write(serialize_edge(from_id, to_id))
