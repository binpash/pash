import argparse
import sys
import socket
import pickle
import traceback
from datetime import datetime
from typing import List, Set, Tuple, Dict, Callable
from uuid import uuid4

sys.path.append("/pash/compiler")

import config
from ir import *
from ast_to_ir import compile_asts
from ir_to_ast import to_shell
from util import *
from dspash.hdfs_utils import HDFSFileConfig

from definitions.ir.aggregator_node import *

from definitions.ir.nodes.eager import *
from definitions.ir.nodes.pash_split import *

import definitions.ir.nodes.r_merge as r_merge
import definitions.ir.nodes.r_split as r_split
import definitions.ir.nodes.r_unwrap as r_unwrap
import definitions.ir.nodes.dgsh_tee as dgsh_tee
import definitions.ir.nodes.remote_pipe as remote_pipe
import shlex
import subprocess
import pash_compiler
from collections import deque, defaultdict

HOST = socket.gethostbyname(socket.gethostname())
NEXT_PORT = 58000
DISCOVERY_PORT = 50052


def read_graph(filename):
    with open(filename, "rb") as ir_file:
        ir, shell_vars = pickle.load(ir_file)
    return ir, shell_vars


def save_configs(graph: IR, dfs_configs_paths: Dict[HDFSFileConfig, str]):
    for edge in graph.all_fids():
        if isinstance(edge.get_resource(), DFSSplitResource):
            resource: DFSSplitResource = edge.get_resource()
            config: HDFSFileConfig = resource.config
            if config not in dfs_configs_paths:
                config_path = ptempfile()
                with open(config_path, "w") as f:
                    f.write(config)
                dfs_configs_paths[config] = config_path
            else:
                config_path = dfs_configs_paths[config]

            resource.set_config_path(config_path)


def to_shell_file(graph: IR, args) -> str:
    filename = ptempfile()

    dirs = set()
    for edge in graph.all_fids():
        directory = os.path.join(config.PASH_TMP_PREFIX, edge.prefix)
        dirs.add(directory)

    for directory in dirs:
        os.makedirs(directory, exist_ok=True)

    if not args.no_eager:
        graph = pash_compiler.add_eager_nodes(graph)

    script = to_shell(graph, args)
    with open(filename, "w") as f:
        f.write(script)
    return filename


def split_ir(graph: IR) -> Tuple[List[IR], Dict[int, IR]]:
    """ Takes an optimized IR and splits it subgraphs. Every subgraph
    is a continues section between a splitter and a merger.
    
    Example: given the following optimized IR
                      - tr -- grep -
                    /                \
    cat - uniq -split                    cat - wc
                    \                /
                      - tr -- grep -
    
    The function returns the following list of IRs:
    [cat - uniq - split, tr - grep, tr - grep, cat - wc]

    Args:
        graph: an IR optimized by the pash compiler
    Returns:
        subgraphs: a list of IRs representing sections of the original graph
        input_fifo_map: a mapping from input edge id to the subgraph which
            have that input edge. This is used later to directly traverse
            from output edge of one subgraph to corrosponding input subgraph

    """
    source_node_ids = graph.source_nodes()
    input_fifo_map = defaultdict(list)

    subgraphs = []
    queue = deque([(source, IR({}, {})) for source in source_node_ids])

    # Graph is a DAG so we need to keep track of traversed edges
    visited_edges = set(graph.all_input_fids())
    visited_nodes = set()

    while queue:
        old_node_id, subgraph = queue.popleft()
        input_fids = graph.get_node_input_fids(old_node_id)
        output_fids = graph.get_node_output_fids(old_node_id)

        if any(map(lambda fid: fid not in visited_edges, input_fids)):
            if subgraph.source_nodes():
                subgraphs.append(subgraph)
            continue

        # Second condition makes sure we don't add empty graphs
        if len(input_fids) > 1 and subgraph.source_nodes():  # merger node
            if subgraph not in subgraphs:
                subgraphs.append(subgraph)
            subgraph = IR({}, {})

        if old_node_id in visited_nodes:
            continue
        else:
            visited_nodes.add(old_node_id)

        node = graph.get_node(old_node_id).copy()
        node_id = node.get_id()

        for idx, input_fid in enumerate(input_fids):
            input_edge_id = None
            # If subgraph is empty and edge isn't ephemeral the edge needs to be added
            if not input_fid.get_ident() in subgraph.edges:
                new_fid = input_fid
                subgraph.add_to_edge(new_fid, node_id)
                input_edge_id = new_fid.get_ident()
            else:
                input_edge_id = input_fid.get_ident()
                subgraph.set_edge_to(input_edge_id, node_id)
            # keep track
            input_fifo_map[input_edge_id].append(subgraph)

        # Add edges coming out of the node
        for output_fid in output_fids:
            subgraph.add_from_edge(node_id, output_fid)
            visited_edges.add(output_fid)

        # Add edges coming into the node
        for input_fid in input_fids:
            if input_fid.get_ident() not in subgraph.edges:
                subgraph.add_to_edge(input_fid, node_id)

        # Add the node
        subgraph.add_node(node)

        next_ids = graph.get_next_nodes(old_node_id)
        if len(next_ids) == 1:
            queue.append((next_ids[0], subgraph))
        else:
            subgraphs.append(subgraph)
            for next_id in next_ids:
                queue.append((next_id, IR({}, {})))

    # print(list(map(lambda k : k.all_fids(), graphs)))
    return subgraphs, input_fifo_map


def add_stdout_fid(graph: IR, file_id_gen: FileIdGen) -> FileId:
    stdout = file_id_gen.next_file_id()
    stdout.set_resource(FileDescriptorResource(("fd", 1)))
    graph.add_edge(stdout)
    return stdout


def assign_workers_to_subgraphs(
    subgraphs: List[IR],
    file_id_gen: FileIdGen,
    input_fifo_map: Dict[int, IR],
    get_worker: Callable,
) -> (IR, Tuple):
    """Takes a list of subgraphs and assigns a worker to each subgraph and augment
    the subgraphs with the necessary remote read/write nodes for data movement
    between workers. This function also produces graph that should run in
    the original shell in which pash was executed. This graph contains
    remote read/write nodes for stdin/stdout, named pipes, and files.

    Args:
        subgraphs: list of sub sections of an optimized IR (returned from split_ir)
        file_id_gen: file id generator of the original ir
        input_fifo_map: mapping from input idge id to subgraph (returned from split_ir)
        get_worker: a callback for getting a worker from worker manager
    Returns:
        main_graph: the graph to execute on main shell
        worker_subgraph_pairs: A list of pairs representing which worker
            each subgraph should be executed on.
    """
    # The graph to execute in the main pash_compiler
    main_graph = IR({}, {})
    worker_subgraph_pairs = []

    # Replace output edges and corrosponding input edges with remote read/write
    for subgraph in subgraphs:
        subgraph_critical_fids = list(
            filter(lambda fid: fid.has_remote_file_resource(), subgraph.all_fids())
        )
        worker = get_worker(subgraph_critical_fids)
        worker._running_processes += 1
        worker_subgraph_pairs.append((worker, subgraph))
        sink_nodes = subgraph.sink_nodes()
        assert len(sink_nodes) == 1

        for out_edge in subgraph.get_node_output_fids(sink_nodes[0]):
            stdout = add_stdout_fid(subgraph, file_id_gen)
            out_edge_id = out_edge.get_ident()
            # Replace the old edge with an ephemeral edge in case it isn't and
            # to avoid modifying the edge in case it's used in some other subgraph
            ephemeral_edge = file_id_gen.next_ephemeral_file_id()
            subgraph.replace_edge(out_edge_id, ephemeral_edge)
            edge_uid = uuid4()
            # Add remote-write node at the end of the subgraph
            remote_write = remote_pipe.make_remote_pipe(
                [ephemeral_edge.get_ident()],
                [stdout.get_ident()],
                worker.host(),
                DISCOVERY_PORT,
                False,
                edge_uid,
            )
            subgraph.add_node(remote_write)

            # Copy the old output edge resource
            new_edge = file_id_gen.next_file_id()
            new_edge.set_resource(out_edge.get_resource())
            # Get the subgraph which "edge" writes to
            if out_edge_id in input_fifo_map and out_edge.is_ephemeral():
                # Copy the old output edge resource
                matching_subgraph = input_fifo_map[out_edge_id][0]
                matching_subgraph.replace_edge(out_edge.get_ident(), new_edge)
            else:
                matching_subgraph = main_graph
                matching_subgraph.add_edge(new_edge)

            remote_read = remote_pipe.make_remote_pipe(
                [],
                [new_edge.get_ident()],
                worker.host(),
                DISCOVERY_PORT,
                True,
                edge_uid,
            )
            matching_subgraph.add_node(remote_read)

    # Replace non ephemeral input edges with remote read/write
    for subgraph in subgraphs:
        source_nodes = subgraph.source_nodes()
        for source in source_nodes:
            for in_edge in subgraph.get_node_input_fids(source):
                if (
                    in_edge.has_file_resource()
                    or in_edge.has_file_descriptor_resource()
                ):
                    # setup
                    stdout = add_stdout_fid(main_graph, file_id_gen)

                    # Copy the old input edge resource
                    new_edge = file_id_gen.next_file_id()
                    new_edge.set_resource(in_edge.get_resource())
                    main_graph.add_edge(new_edge)

                    # Add remote write to main subgraph
                    edge_uid = uuid4()
                    remote_write = remote_pipe.make_remote_pipe(
                        [new_edge.get_ident()],
                        [stdout.get_ident()],
                        HOST,
                        DISCOVERY_PORT,
                        False,
                        edge_uid,
                    )
                    main_graph.add_node(remote_write)

                    # Add remote read to current subgraph
                    ephemeral_edge = file_id_gen.next_ephemeral_file_id()
                    subgraph.replace_edge(in_edge.get_ident(), ephemeral_edge)

                    remote_read = remote_pipe.make_remote_pipe(
                        [],
                        [ephemeral_edge.get_ident()],
                        HOST,
                        DISCOVERY_PORT,
                        True,
                        edge_uid,
                    )
                    subgraph.add_node(remote_read)
                else:
                    # sometimes a command can have both a file resource and an ephemeral resources (example: spell oneliner)
                    continue

    return main_graph, worker_subgraph_pairs


def prepare_graph_for_remote_exec(filename: str, get_worker: Callable):
    """
    Reads the complete ir from filename and splits it
    into subgraphs where ony the first subgraph represent a continues
    segment (merger segment or branched segment) in the graph.
    Note: All subgraphs(except first one) read and write from remote pipes.
        However, we had to add a fake stdout to avoid some problems when converting to shell code.

    Returns:
        worker_graph_pairs: List of (worker, subgraph)
        shell_vars: shell variables
        main_graph: The ir we need to execute on the main shell.
            This graph contains edges to correctly redirect the following to remote workers
            - special pipes (stdin/stdout)
            - named pipes reading and writing
            - files reading and writing
    """
    ir, shell_vars = read_graph(filename)
    file_id_gen = ir.get_file_id_gen()
    subgraphs, mapping = split_ir(ir)
    main_graph, worker_graph_pairs = assign_workers_to_subgraphs(
        subgraphs, file_id_gen, mapping, get_worker
    )
    return worker_graph_pairs, shell_vars, main_graph
