import argparse
import sys
import socket
import pickle
import traceback
import collections
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

from definitions.ir.nodes.eager import *
from definitions.ir.nodes.pash_split import *

import definitions.ir.nodes.r_merge as r_merge
import definitions.ir.nodes.r_split as r_split
import definitions.ir.nodes.r_unwrap as r_unwrap
import definitions.ir.nodes.dgsh_tee as dgsh_tee
import definitions.ir.nodes.remote_pipe as remote_pipe
import shlex
import subprocess
from collections import deque, defaultdict
import stat, os
import pash_compiler

HOST = socket.gethostbyname(socket.gethostname())
NEXT_PORT = 58000
DISCOVERY_PORT = 50052


def read_graph(filename):
    with open(filename, "rb") as ir_file:
        ir, shell_vars = pickle.load(ir_file)
    return ir, shell_vars

def save_configs(graph:IR, dfs_configs_paths: Dict[HDFSFileConfig, str]):
    for edge in graph.all_fids():
        if isinstance(edge.get_resource(), DFSSplitResource):
            resource : DFSSplitResource = edge.get_resource()
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
        # Set DFGNode next id to not clash with already existing ids
        # TODO: ideally we should get the next_id from the graph object
        #   to avoid conflicts across parallel processes
        DFGNode.next_id = max(DFGNode.next_id , max(graph.nodes.keys()) + 1)
        graph = pash_compiler.add_eager_nodes(graph)
    script = to_shell(graph, args)
    with open(filename, "w") as f:
        f.write(script)
    return filename

def optimize_named_fifos(graph):
    """
    Replaces named fifos with ephemeral fifos when we
    know both the read and write ends of the named fifo
    """
    named_fifos = set()
    for fid in graph.all_fids():
        path = str(fid.get_resource())
        try:
            if fid.has_file_resource() and stat.S_ISFIFO(os.stat(path).st_mode):
                _, from_node, to_node = graph.edges[fid.ident]
                if from_node != None and to_node != None:
                    named_fifos.add(fid)
        except:
            pass
    for named_fifo in named_fifos:
        named_fifo.make_ephemeral()

    return graph

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
    
    graph = optimize_named_fifos(graph)
    
    subgraphs = set()
    queue = deque([(source, IR({}, {})) for source in source_node_ids])

    # Set next_graph_policy to the wanted policy
    combine_after_merge_policy = lambda combining_subgraph=queue[0][1]: combining_subgraph
    next_graph_policy = combine_after_merge_policy
    # Comment the above and uncomment below to change graph splitting policy
    # new_graph_policy: lambda : IR({}, {})
    # next_graph_policy = new_graph_policy

    # Graph is a DAG so we need to keep track of traversed edges
    visited_edges = set(graph.all_input_fids())
    visited_nodes = set()

    while queue:
        old_node_id, subgraph = queue.popleft()
        input_fids = graph.get_node_input_fids(old_node_id)
        output_fids = graph.get_node_output_fids(old_node_id)

        if(any(map(lambda fid:fid not in visited_edges, input_fids))):
            if subgraph.source_nodes():
                subgraphs.add(subgraph)
            continue
        
        # Second condition makes sure we don't add empty graphs
        if len(input_fids) > 1 and subgraph.source_nodes(): # merger node
            if subgraph not in subgraphs:
                subgraphs.add(subgraph)
            subgraph = next_graph_policy()

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
            subgraphs.add(subgraph)
            for next_id in next_ids:
                queue.append((next_id, next_graph_policy()))

    # for graph in subgraphs:
    #     print(to_shell(graph, config.pash_args), file=sys.stderr)
     
    return subgraphs, input_fifo_map

def add_stdout_fid(graph : IR, file_id_gen: FileIdGen) -> FileId:
    stdout = file_id_gen.next_file_id()
    stdout.set_resource(FileDescriptorResource(('fd', 1)))
    graph.add_edge(stdout)
    return stdout


def create_remote_pipe_from_output_edge(from_subgraph: IR, to_subgraph: IR, edge: FileId, host, port, file_id_gen: FileIdGen):
    stdout = add_stdout_fid(from_subgraph, file_id_gen)
    edge_id = edge.get_ident()

    # Replace the old edge with an ephemeral edge in case it isn't and
    # to avoid modifying the edge in case it's used in some other subgraph
    ephemeral_edge = file_id_gen.next_ephemeral_file_id()
    from_subgraph.replace_edge(edge_id, ephemeral_edge)
    edge_uid = uuid4()
    # Add remote-write node at the end of the subgraph
    remote_write = remote_pipe.make_remote_pipe([ephemeral_edge.get_ident()], [stdout.get_ident()], host, port, False, edge_uid)
    from_subgraph.add_node(remote_write)
    
    # Copy the old output edge resource
    new_edge = file_id_gen.next_file_id()
    new_edge.set_resource(edge.get_resource())
    # Get the subgraph which "edge" writes to
    if edge_id in to_subgraph.edges:
        # Replace the old output edge resource
        to_subgraph.replace_edge(edge_id, new_edge)
    else:
        to_subgraph.add_edge(new_edge)

    remote_read = remote_pipe.make_remote_pipe([], [new_edge.get_ident()], host, port, True, edge_uid)
    to_subgraph.add_node(remote_read)

def create_remote_pipe_from_input_edge(from_subgraph: IR, to_subgraph: IR, edge: FileId, host, port, file_id_gen: FileIdGen):
    stdout = add_stdout_fid(from_subgraph, file_id_gen)

    # Copy the old input edge resource
    new_edge = file_id_gen.next_file_id()
    new_edge.set_resource(edge.get_resource())
    from_subgraph.add_edge(new_edge)

    # Add remote write to main subgraph
    edge_uid = uuid4()
    remote_write = remote_pipe.make_remote_pipe([new_edge.get_ident()], [stdout.get_ident()], HOST, DISCOVERY_PORT, False, edge_uid)
    from_subgraph.add_node(remote_write)

    # Add remote read to current subgraph
    ephemeral_edge = file_id_gen.next_ephemeral_file_id()
    old_edge_id = edge.get_ident()

    to_subgraph.replace_edge(old_edge_id, ephemeral_edge)

    remote_read = remote_pipe.make_remote_pipe([], [ephemeral_edge.get_ident()], HOST, DISCOVERY_PORT, True, edge_uid)
    to_subgraph.add_node(remote_read)

def assign_workers_to_subgraphs(subgraphs:List[IR], file_id_gen: FileIdGen, input_fifo_map:Dict[int, IR], get_worker: Callable) -> Tuple[IR, List]:
    """ Takes a list of subgraphs and assigns a worker to each subgraph and augment
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
        subgraph_critical_fids = list(filter(lambda fid: fid.has_remote_file_resource(), subgraph.all_fids()))
        worker = get_worker(subgraph_critical_fids)
        worker._running_processes += 1
        worker_subgraph_pairs.append([worker, subgraph])
        sink_nodes = subgraph.sink_nodes()
        
        for sink_node in sink_nodes:
            for out_edge in subgraph.get_node_output_fids(sink_node):
                out_edge_id = out_edge.get_ident()
                if out_edge_id in input_fifo_map and out_edge.is_ephemeral():
                    # Copy the old output edge resource
                    matching_subgraph = input_fifo_map[out_edge_id][0]
                else:
                    matching_subgraph = main_graph

                if matching_subgraph == subgraph:
                    continue
                
                # In case this fid from and to nodes are in the same graph
                # we need to both read/write the fid from/to the main machine
                # TODO: we might want to optimize this by offloading everything to the
                # main subgraph
                to_node = subgraph.get_edge_to(out_edge_id)
                if to_node:
                    # remove previous connection
                    subgraph.set_edge_to(out_edge_id, None)

                    # Add new connectiong 
                    new_output_edge = file_id_gen.next_file_id()
                    new_output_edge.set_resource(out_edge.get_resource())
                    subgraph.add_to_edge(new_output_edge, to_node)
                    subgraph.get_node(to_node).replace_edge(out_edge_id, new_output_edge.get_ident())

                    # Add remote pipe to write from subgraph to main graph
                    create_remote_pipe_from_input_edge(from_subgraph=matching_subgraph, to_subgraph=subgraph, edge=new_output_edge,  host=HOST, port=DISCOVERY_PORT, file_id_gen=file_id_gen)


                create_remote_pipe_from_output_edge(from_subgraph=subgraph, to_subgraph=matching_subgraph, edge=out_edge, host=worker.host(), port=DISCOVERY_PORT, file_id_gen=file_id_gen)

    # Replace non ephemeral input edges with remote read/write
    for worker, subgraph in worker_subgraph_pairs:
        nodes = list(subgraph.nodes.keys())
        for source in nodes:
            for in_edge in subgraph.get_node_input_fids(source):
                if in_edge.has_file_resource() or in_edge.has_file_descriptor_resource():
                    old_edge_id = in_edge.get_ident()

                    # In case this fid from and to nodes are in the same graph
                    # we need to both read/write the fid from/to the main machine
                    # TODO: we might want to optimize this by offloading everything to the
                    # main subgraph
                    from_node = subgraph.get_edge_from(old_edge_id)
                    if from_node:
                        # remove previous connection
                        subgraph.set_edge_from(old_edge_id, None)

                        # Add new connectiong 
                        new_output_edge = file_id_gen.next_file_id()
                        new_output_edge.set_resource(in_edge.get_resource())
                        subgraph.add_from_edge(from_node, new_output_edge)
                        subgraph.get_node(from_node).replace_edge(old_edge_id, new_output_edge.get_ident())

                        # Add remote pipe to write from subgraph to main graph
                        create_remote_pipe_from_output_edge(from_subgraph=subgraph, to_subgraph=main_graph, edge=new_output_edge,  host=worker.host(), port=DISCOVERY_PORT, file_id_gen=file_id_gen)

                    create_remote_pipe_from_input_edge(from_subgraph=main_graph, to_subgraph=subgraph, edge=in_edge,  host=HOST, port=DISCOVERY_PORT, file_id_gen=file_id_gen)
                else:
                    # sometimes a command can have both a file resource and an ephemeral resources (example: spell oneliner)
                    continue

    # for worker, graph in worker_subgraph_pairs:
    #     print(to_shell(graph, config.pash_args), file=sys.stderr)

    return main_graph, worker_subgraph_pairs

def add_debug_flags(graph: IR):
    for node in graph.nodes.values():
        if isinstance(node, remote_pipe.RemotePipe):
            node.add_debug_flag()

def prepare_graph_for_remote_exec(filename:str, get_worker:Callable):
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
    main_graph, worker_graph_pairs = assign_workers_to_subgraphs(subgraphs, file_id_gen, mapping, get_worker)
    return worker_graph_pairs, shell_vars, main_graph, file_id_gen, mapping

def get_best_worker_for_subgraph(subgraph:IR, get_worker:Callable):
    """Get the best worker for a subgraph
    
    subgraph: the ir we want to match with the best worker
    get_worker: a callback for getting a worker from worker manager
    Return: the best worker
    """    
    subgraph_critical_fids = list(filter(lambda fid: fid.has_remote_file_resource(), subgraph.all_fids()))
    worker = get_worker(subgraph_critical_fids)
    worker._running_processes += 1
    return worker

def get_neighbor_remote_reader_pipes(subgraphs: [IR], remote_pipe_id, exclude_subgraph: IR) -> [remote_pipe.RemotePipe]:
    """Get all (subgraph, node) pairs such that the node has established communication
        with the crashed_worker. This includes remote_pipes on the crashed_worker too

    Assumption: the subgraph has already gone through the initial processing in prepare_graph_for_remote_exec()
                therefore, all remote_writer_nodes are sink_nodes 
                and all remote_reader_nodes are source_nodes 
                (not true the other way - source_nodes can also be DFSSplitReader nodes)

    subgraphs: list of IRs after split
    remote_pipe_id: id of the remote pipe we want to find neighbors for
    Return: set of nodes whose addresses need to be updated
    """
    update_remote_pipe_candidates = []
    for subgraph in subgraphs:
        if subgraph != exclude_subgraph:
            for source_id in subgraph.source_nodes():
                source_node = subgraph.get_node(source_id)
                # Sanity check
                if isinstance(source_node, remote_pipe.RemotePipe):
                    assert(source_node.is_remote_read())
                    # If the current remote_pipe had communicated with remote_pipe_id before,
                    # append it to the set
                    if source_node.get_uuid() == remote_pipe_id:
                        update_remote_pipe_candidates.append(source_node)
    return update_remote_pipe_candidates

def get_remote_pipe_update_candidates(subgraphs: [IR], subgraph: IR, crashed_worker):
    """Update remote nodes (remote writer/reader)'s host and port fields
       with the new neighbor's host

    Assumption: the subgraph has already gone through the initial processing in prepare_graph_for_remote_exec()
                therefore, all remote_writer_nodes are sink_nodes 
                and all remote_reader_nodes are source_nodes 
                (not true the other way - source_nodes can also be DFSSplitReader nodes)

    subgraphs: list of IRs after split
    subgraph: the target subgraph (belonged to the crashed worker) whose remote_pipe addresses [rp] we need to update
              we also need to update remote_pipe addresses for any neighbor remote_pipes for each of the [rp] 
    crashed_worker: the original crashed worker (type: WorkerConnection)
    replacement_worker: the replacement worker (type: WorkerConnection) whose host should be used to update the nodes
    Return: {update_candidate:replacement_worker} where update_candidate is any remote pipe candidate for update that belong to subgraph
    """
    # Find all remote_pipes for every remote_pipe in subgraph that needs to be updated
    # remote_pipes whose addr is not the crashed worker don't need to be updated 
    # (their addresses will be updated to be the same replacement_worker)
    update_candidates = []
    for boundary_node_id in subgraph.sink_nodes():
        boundary_node = subgraph.get_node(boundary_node_id)
        if isinstance(boundary_node, remote_pipe.RemotePipe) and boundary_node.get_host() == crashed_worker.host():
            assert(boundary_node.is_remote_read() == False)
            update_candidates.append(boundary_node)
            neighbor_reader_rp = get_neighbor_remote_reader_pipes(subgraphs, boundary_node.get_uuid(), subgraph)
            update_candidates = update_candidates + neighbor_reader_rp

    return update_candidates


def update_remote_pipe_addr(update_node, new_host, new_port=DISCOVERY_PORT):
    assert(isinstance(update_node, remote_pipe.RemotePipe))
    update_node.set_addr(host_ip=new_host, port=new_port)

def update_subgraphs(original_to_updated_subgraphs, update_candidates_replacements_map):
    for original_subgraph in original_to_updated_subgraphs.keys():
        new_subgraph = copy.deepcopy(original_subgraph)
        for boundary_node_id in new_subgraph.source_nodes() + new_subgraph.sink_nodes():
            boundary_node = new_subgraph.get_node(boundary_node_id)
            if isinstance(boundary_node, remote_pipe.RemotePipe):
                if boundary_node in update_candidates_replacements_map:
                    update_remote_pipe_addr(boundary_node, update_candidates_replacements_map[boundary_node])
        
        original_to_updated_subgraphs[original_subgraph] = new_subgraph

def get_worker_subgraph_map(worker_subgraph_pairs):
    worker_subgraph_map = collections.defaultdict(list)
    for worker, subgraph in worker_subgraph_pairs:
        worker_subgraph_map[worker].append(subgraph)
    return worker_subgraph_map