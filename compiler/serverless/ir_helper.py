# splitter, merger: change key to key_v0, fifo in the pashlib to fifo_v0, for each pashlib, know the downstream script id, add command
# lambda: change key to key_v${version}

import argparse
from copy import deepcopy
from fileinput import filename
from typing import Dict, List, Tuple
import sys
import os
import json
from uuid import uuid4

import boto3
import time
from contextlib import contextmanager

sys.path.append(os.path.join(os.getenv("PASH_TOP"), "compiler"))
from compiler.serverless.graph_print_helper import pretty_print_subgraphs
import definitions.ir.nodes.serverless_remote_pipe as serverless_remote_pipe
import definitions.ir.nodes.serverless_lambda_invoke as serverless_lambda_invoke
from definitions.ir.nodes.r_wrap import RWrap
from definitions.ir.nodes.r_split import RSplit
from definitions.ir.nodes.r_merge import RMerge
from definitions.ir.nodes.cat import make_cat_node
from dspash.ir_helper import split_ir
from ir_to_ast import to_shell
from ir import *
import config
import pash_compiler
from pash_annotations.datatypes.AccessKind import make_stream_output, make_stream_input
from pash_annotations.datatypes.BasicDatatypes import Operand
from pash_annotations.datatypes.CommandInvocationWithIOVars import CommandInvocationWithIOVars
from definitions.ir.dfg_node import DFGNode
from definitions.ir.arg import Arg


# ============================================================================
# Timing instrumentation utilities
# ============================================================================

# Global timing store
_timing_data = {}

@contextmanager
def time_block(label):
    """Context manager to time a code block and store results."""
    start = time.time()
    try:
        yield
    finally:
        elapsed = time.time() - start
        _timing_data[label] = elapsed
        print(f"[TIMING] {label}: {elapsed:.3f}s")

def get_timing(label):
    """Retrieve a timing measurement."""
    return _timing_data.get(label, 0.0)

def print_timing_summary():
    """Print all timing measurements."""
    if not _timing_data:
        return
    print("\n" + "="*60)
    print("TIMING SUMMARY")
    print("="*60)
    for label, elapsed in _timing_data.items():
        print(f"  {label:<40} {elapsed:>8.3f}s")
    print("="*60 + "\n")


# Debugpy setup - connect before main logic runs
debug = False
import debugpy
if debug:
    debugpy.listen(5000)
    print("🐛 Debugger listening on port 5684. Attach VSCode now!")
    debugpy.wait_for_client()
    print("✅ Debugger attached! Continuing execution...")

graph = False

def add_stdout_fid(graph : IR, file_id_gen: FileIdGen) -> FileId:
    stdout = file_id_gen.next_file_id()
    stdout.set_resource(FileDescriptorResource(('fd', 1)))
    graph.add_edge(stdout)
    return stdout


def adjust_lambda_incoming_edges(
    first_subgraph: IR,
    subgraphs: List[IR],
    file_id_gen: FileIdGen,
    s3_file_path: str,
    rsplit_output_ids: List[int],
    ec2_file_in_edge: FileId
) -> None:
    """
    Adjust incoming edges for lambda subgraphs after removing an initial ServerlessRemotePipe
    used in the S3 -> cat -> split direct streaming optimization.

    Replaces rsplit output edges with ServerlessRemotePipe nodes that read directly from S3.
    All downstream subgraphs read from the same S3 file (not partitioned).

    Args:
        first_subgraph: the (modified) first subgraph after removing the remote-read node.
        subgraphs: the list of all subgraphs produced by split_ir (may be mutated).
        file_id_gen: generator to allocate new FileId / ephemeral ids.
        s3_file_path: S3 file path to read from (e.g., "input.txt" or "s3://bucket/file.txt")
        rsplit_output_ids: List of edge IDs that were outputs from the rsplit node
    """


    #
    """
    
    input/1M.txt -> cat -> split            fifo1 -> lambda1 (sort) -> ...
                                             fifo -> lambda2 (sort) -> ...
                                            
                                            input/1M.txt -> sort
                                             
                            
    

    
    """
    # Convert to set as we only need to do this once per unique edge/lambda
    rsplit_output_id_set = set(rsplit_output_ids)
    total_lambdas = 0

    # Iterate through all downstream subgraphs
    for subgraph in subgraphs:
        source_nodes = subgraph.source_nodes()

        for source_node_id in source_nodes:
            source_node = subgraph.get_node(source_node_id)

            # Get all input edges of this source node
            input_fids = subgraph.get_node_input_fids(source_node_id)

            for in_edge in input_fids:
                in_edge_id = in_edge.get_ident()

                # Is this edge one of the rsplit outputs?
                if in_edge_id in rsplit_output_id_set:
                    # YES - Replace with in edge from ec2
                    subgraph.replace_edge(in_edge_id, ec2_file_in_edge)
                    total_lambdas += 1


                    # next step: replace with a byte range command 


                    #TODO. replace incoming edges of lambdas with the in_edge

                    #ephemeral_edge = file_id_gen.next_ephemeral_file_id()

                    # 2. Replace old rsplit edge with new edge
                    #subgraph.replace_edge(in_edge_id, ephemeral_edge)

                    # 3. Add S3 read node
                    # #remote_read = serverless_remote_pipe.make_serverless_remote_pipe(
                    #     local_fifo_id=ephemeral_edge.get_ident(),
                    #     is_remote_read=True,
                    #     remote_key=s3_file_path,
                    #     output_edge=None,
                    #     is_tcp=False
                    # )
                    #subgraph.add_node(remote_read)

                    # 4. Add eager/dgsh-tee after S3 read (THIS IS THE MISSING PIECE!)
                    #pash_compiler.add_eager(ephemeral_edge.get_ident(), subgraph, file_id_gen)


                    """
                    # 1. Create ephemeral edge
                    file_id_gen.next_ephemeral_file_id()
                    ephemeral_edge = file_id_gen.next_ephemeral_file_id() #todo. is this not unique 
                    # 2. Replace old edge (adds ephemeral_edge to graph, connects TO source_node)
                    subgraph.replace_edge(in_edge_id, ephemeral_edge)
                    #subgraph.add_to_edge(ephemeral_edge, source_node_id)

                    # 3. Create remote_read node
                    remote_read = serverless_remote_pipe.make_serverless_remote_pipe(
                        local_fifo_id=ephemeral_edge.get_ident(),  # <-  Sets as implicit output
                        is_remote_read=True,
                        remote_key=s3_file_path,
                        output_edge=None,
                        is_tcp=False
                    )

                    # 4. Add node (auto-connects remote_read → ephemeral_edge via set_edge_from)
                    subgraph.add_node(remote_read)
                    #then we need to connect it to the source node 
                    # we need to create a new edge for this 
                    ephemeral_edge2 = file_id_gen.next_ephemeral_file_id()

                    subgraph.add_from_edge(remote_read.get_id(), ephemeral_edge2)
                    subgraph.set_edge_to(ephemeral_edge2.get_ident(), source_node_id)
                    """
                    #TODO. make this better or use the implicit ouput so you don't have to create a new edge like this

                    """
                    # 1. Create new ephemeral edge
                    ephemeral_edge = file_id_gen.next_ephemeral_file_id()

                    # 2. Add ephemeral edge to graph
                    subgraph.add_edge(ephemeral_edge)

                    # 3. Create remote_read node
                    remote_read = serverless_remote_pipe.make_serverless_remote_pipe(
                        local_fifo_id=ephemeral_edge.get_ident(),
                        is_remote_read=True,
                        remote_key=s3_file_path,
                        output_edge=None,
                        is_tcp=False
                    )

                    # 4. Add remote_read node (auto-connects its declared edges)
                    subgraph.add_node(remote_read)

                    # 5. Connect remote_read → ephemeral_edge
                    subgraph.add_from_edge(remote_read.get_id(), ephemeral_edge)

                    # 6. Update source_node to use ephemeral_edge as input (targeted, not global)
                    source_node = subgraph.get_node(source_node_id)
                    source_node.replace_edge(in_edge_id, ephemeral_edge.get_ident())

                    # 7. Connect ephemeral_edge → source_node
                    subgraph.set_edge_to(ephemeral_edge.get_ident(), source_node_id)

                    # 8. Disconnect old edge from source_node
                    subgraph.set_edge_to(in_edge_id, None)
                    """
    return total_lambdas


# ============================================================================
# Byte-range optimization helpers (for S3 direct lambda streaming)
# ============================================================================

def unwrap_rwrap_node(rwrap_node: RWrap) -> DFGNode:
    """
    Extract the inner command from an RWrap node and create a plain DFGNode
    that runs it directly without the r_wrap wrapper.

    RWrap structure: r_wrap bash -c 'command'
    Unwrapped:       bash -c 'command'
    """
    cmd_inv = rwrap_node.cmd_invocation_with_io_vars

    # operand_list = ["bash -c", "'command args...'"]
    # We want to run: bash -c 'command args...'
    operand_list = cmd_inv.operand_list

    # Extract the inner command string (operand_list[1])
    inner_cmd = operand_list[1] if len(operand_list) > 1 else operand_list[0]

    input_id = cmd_inv.implicit_use_of_streaming_input
    output_id = cmd_inv.implicit_use_of_streaming_output

    access_map = {
        input_id: make_stream_input(),
        output_id: make_stream_output()
    }

    # Create bash -c command
    new_cmd_inv = CommandInvocationWithIOVars(
        cmd_name="bash",
        flag_option_list=[],
        operand_list=[Operand(Arg.string_to_arg("-c")), inner_cmd],
        implicit_use_of_streaming_input=input_id,
        implicit_use_of_streaming_output=output_id,
        access_map=access_map
    )

    return DFGNode(new_cmd_inv,
                   com_redirs=rwrap_node.com_redirs,
                   com_assignments=rwrap_node.com_assignments)


def unwrap_rwraps_in_subgraph(subgraph: IR) -> int:
    """
    Find all RWrap nodes in a subgraph and replace them with unwrapped
    versions that run the inner command directly.

    Returns: Number of RWrap nodes unwrapped
    """
    unwrapped_count = 0

    # Collect RWrap node IDs first (avoid modifying dict during iteration)
    rwrap_ids = [nid for nid, node in subgraph.nodes.items() if isinstance(node, RWrap)]

    for rwrap_id in rwrap_ids:
        rwrap_node = subgraph.get_node(rwrap_id)

        # Create unwrapped node with same edges
        unwrapped_node = unwrap_rwrap_node(rwrap_node)

        # Remove old RWrap node (disconnects edges)
        subgraph.remove_node(rwrap_id)

        # Add new unwrapped node (reconnects edges)
        subgraph.add_node(unwrapped_node)

        unwrapped_count += 1

    return unwrapped_count


def replace_rmerge_with_cat(subgraph: IR) -> bool:
    """
    Replace RMerge nodes (runtime/r_merge) with cat nodes.

    Returns: True if any replacement was made
    """
    replaced = False

    # Collect RMerge node IDs
    rmerge_ids = [nid for nid, node in subgraph.nodes.items() if isinstance(node, RMerge)]

    for rmerge_id in rmerge_ids:
        rmerge_node = subgraph.get_node(rmerge_id)

        input_ids = rmerge_node.get_input_list()
        output_ids = rmerge_node.get_output_list()

        if len(output_ids) != 1:
            continue  # Unexpected structure

        output_id = output_ids[0]

        # Create cat node with same inputs and output
        cat_node = make_cat_node(input_ids, output_id)

        # Remove old RMerge
        subgraph.remove_node(rmerge_id)

        # Add cat node (this reconnects the edges)
        subgraph.add_node(cat_node)

        replaced = True

    return replaced


def is_merger_subgraph(subgraph: IR) -> bool:
    """
    Check if a subgraph is a merger/EC2 node by looking for RMerge nodes.
    """
    for node in subgraph.nodes.values():
        if isinstance(node, RMerge):
            return True
    return False


def get_downstream_subgraphs(subgraph: IR, input_fifo_map: Dict[int, Tuple]) -> List[IR]:
    """
    Find all subgraphs that consume outputs from the given subgraph.

    Args:
        subgraph: The source subgraph
        input_fifo_map: Mapping from output edge id to (consuming_subgraph, ...)

    Returns:
        List of downstream subgraphs
    """
    downstream = []
    sink_nodes = subgraph.sink_nodes()

    for sink_id in sink_nodes:
        out_edges = subgraph.get_node_output_fids(sink_id)
        for out_edge in out_edges:
            out_edge_id = out_edge.get_ident()
            if out_edge_id in input_fifo_map:
                downstream_subgraph = input_fifo_map[out_edge_id][0]
                if downstream_subgraph not in downstream:
                    downstream.append(downstream_subgraph)

    return downstream


"""
Apply byte-range optimizations across a set of subgraphs reachable from given worker subgraphs.

    This function performs a breadth-first traversal starting from the subgraphs listed in
    `worker_subgraphs` and applies two byte-range related optimizations to each visited subgraph:
    - unwrap RWrap nodes (collapsing unnecessary wrappers around byte-range values)
    - replace RMerge nodes in merger subgraphs with a concatenation ("cat") operation

    Traversal and behavior:
    - The traversal is BFS and tracks visited subgraphs by their Python id() to avoid revisiting.
    - For each visited subgraph:
        1. Call unwrap_rwraps_in_subgraph(subgraph). If any RWrap nodes are unwrapped, a message
             is printed summarizing how many were removed.
        2. If is_merger_subgraph(subgraph) returns True, the function attempts to replace any
             RMerge node(s) in that subgraph with a cat operation by calling replace_rmerge_with_cat(subgraph).
             If replacement occurs, a message is printed. After processing a merger subgraph, traversal
             does NOT continue past it (the merger is a stopping point).
        3. If the subgraph is not a merger, downstream subgraphs are discovered via
             get_downstream_subgraphs(subgraph, input_fifo_map) and enqueued for processing (unless already visited).

    Parameters:
    - worker_subgraphs (List[IR]):
            List of IR subgraph objects that are the initial workers fed by rsplit outputs.
            These serve as BFS roots for the optimization pass.
    - all_subgraphs (List[IR]):
            All subgraphs in the pipeline. Note: the current implementation does not use this
            argument directly, but it can be provided for future needs or for callers that
            expect a complete pipeline view.
    - input_fifo_map (Dict[int, Tuple]):
            Mapping used by get_downstream_subgraphs to resolve downstream connectivity. Keys are
            FIFO identifiers (or other integer keys) and values are tuples describing the consumer
            connection. The exact shape is determined by the graph representation and helper utilities.

"""
def apply_byte_range_optimizations(
    
    worker_subgraphs: List[IR],
    all_subgraphs: List[IR],
    input_fifo_map: Dict[int, Tuple]
) -> None:
    """
    Apply byte-range optimizations (unwrap RWrap, replace RMerge) to:
    1. Worker subgraphs (fed by rsplit outputs)
    2. Their consumer subgraphs (may be multiple consecutive workers)
    3. Continue until reaching a merger subgraph (contains RMerge)
    4. Apply to merger, then STOP (don't process downstream of merger)

    Example pipeline:
      rsplit -> worker1 -> worker2 -> worker3 -> merger (RMerge) -> downstream
                ↑          ↑          ↑          ↑
                process    process    process    process & STOP
    """
    processed = set()  # Track processed subgraphs by id
    queue = list(worker_subgraphs)  # BFS queue

    while queue:
        subgraph = queue.pop(0)
        sg_id = id(subgraph)

        if sg_id in processed:
            continue
        processed.add(sg_id)

        # Apply optimizations to this subgraph
        rwraps_unwrapped = unwrap_rwraps_in_subgraph(subgraph)
        if rwraps_unwrapped > 0:
            print(f"[IR Helper] Unwrapped {rwraps_unwrapped} RWrap nodes")

        # Check if this is a merger subgraph
        if is_merger_subgraph(subgraph):
            # Replace RMerge with cat in the merger
            if replace_rmerge_with_cat(subgraph):
                print("[IR Helper] Replaced RMerge with cat in merger subgraph")
            # STOP - don't process downstream of merger
            continue

        # Not a merger - continue to downstream subgraphs
        downstream = get_downstream_subgraphs(subgraph, input_fifo_map)
        for ds in downstream:
            if id(ds) not in processed:
                queue.append(ds)


# cat -> split
def optimize_s3_lambda_direct_streaming(subgraphs:List[IR], input_fifo_map: Dict[int, Tuple] = None):
    if len(subgraphs) == 1: #no point in optimizing here 
        return subgraphs, None, 0
    
    first_subgraph = subgraphs[0]
    file_id_gen = first_subgraph.get_file_id_gen()

    #1. check incoming edge for file resource
    #if true then check if node is cat
    # if true then check if next node is split
    # if true then we can remove and
    # at any point if false return subgraphs

    # byte range 0-500 500-1000
    #worst case have to read entire file
    # expected best case, lambdas can coordinate to read different parts of the file
    #if "file->cat->split":
        #change lambda incoming edges
        #pash posh shark 

    source_nodes = first_subgraph.source_nodes()    # list of ints
    if len(source_nodes) != 1:
        return subgraphs, None, 0
    #TODO. encode assertions to run on narrow execution path 
    for source in source_nodes:
        in_edges = first_subgraph.get_node_input_fids(source)
        if len(in_edges) != 1:
            return subgraphs, None, 0
        
        for in_edge in in_edges:
            # (isinstance(self.resource, FileResource))
            if in_edge.has_file_resource() or in_edge.has_file_descriptor_resource():
                source_node = first_subgraph.get_node(source)
                if hasattr(source_node, 'cmd_invocation_with_io_vars') and source_node.cmd_invocation_with_io_vars.cmd_name == 'cat':
                    #2. find if cat
                    # Extract S3 file path from cat input
                    s3_file_path = in_edge.get_resource().uri

                    next_nodes_ids = first_subgraph.get_next_nodes(source)

                    #           ir obj has no attribute get edge consumers

                    if len(next_nodes_ids) != 1:
                        return subgraphs, None, 0
                    next_node_id = next_nodes_ids[0]
                    next_node = first_subgraph.get_node(next_node_id)

                    if isinstance(next_node, RSplit):
                        # Extract rsplit output edge IDs
                        #can do in a more straightforward way with a helper method in ir
                        rsplit_output_fids = first_subgraph.get_node_output_fids(next_node_id)
                        rsplit_output_ids = [fid.get_ident() for fid in rsplit_output_fids]

                        # Call adjust function with extracted info
                        # TODO. dont need first subgraph anymore
                        total_lambdas = adjust_lambda_incoming_edges(
                            first_subgraph,
                            subgraphs[1:],
                            file_id_gen,
                            s3_file_path,       # S3 file path
                            rsplit_output_ids,  # rsplit output edge IDs,
                            in_edge
                        )

                        # Apply byte-range optimizations ONLY in single-chunk mode
                        chunks_per_lambda_opt = int(os.environ.get('PASH_S3_CHUNKS_PER_LAMBDA', '1'))

                        # Apply byte-range optimizations (unwrap RWrap, replace RMerge with cat)
                        # from workers through to the merger subgraph
                        if total_lambdas > 0 and input_fifo_map is not None:
                            # Worker subgraphs are those that receive rsplit output edges
                            worker_subgraphs = []
                            for rsplit_output_id in rsplit_output_ids:
                                if rsplit_output_id in input_fifo_map:
                                    worker_sg = input_fifo_map[rsplit_output_id][0]
                                    if worker_sg not in worker_subgraphs:
                                        worker_subgraphs.append(worker_sg)

                            # Apply optimizations from workers through to merger
                            if worker_subgraphs:
                                if chunks_per_lambda_opt == 1:
                                    # Single-chunk mode: Apply optimizations (unwrap RWrap, replace r_merge with cat)
                                    apply_byte_range_optimizations(worker_subgraphs, subgraphs, input_fifo_map)
                                    print("[IR Helper] Applied byte-range optimizations (unwrapped RWrap, replaced r_merge)")
                                else:
                                    # Multi-chunk mode: SKIP optimizations - we NEED r_merge for block ordering
                                    print(f"[IR Helper] Skipping byte-range optimizations (chunks_per_lambda={chunks_per_lambda_opt})")
                                    print("[IR Helper] Keeping RWrap nodes and r_merge for multi-chunk block ordering")

                        return subgraphs[1:], in_edge, total_lambdas

    return subgraphs, None, 0




def expand_search_for_newline(s3, bucket, key, start_pos, file_size, max_search=100*1024*1024):
    """
    Exponentially expand search window if newline not found in initial window.
    Handles edge cases like very long lines (>1MB).

    Args:
        s3: boto3 S3 client
        bucket: S3 bucket name
        key: S3 key
        start_pos: Starting byte position to search from
        file_size: Total file size
        max_search: Maximum bytes to search (default 100MB)

    Returns:
        Tuple of (position, found) where:
        - position: Byte position of newline +1, or end of search range if not found
        - found: True if newline was found, False otherwise
    """
    current_pos = start_pos
    chunk_size = 1 * 1024 * 1024  # Start with 1MB

    while current_pos < min(start_pos + max_search, file_size):
        try:
            end_pos = min(current_pos + chunk_size - 1, file_size - 1)
            response = s3.get_object(
                Bucket=bucket,
                Key=key,
                Range=f'bytes={current_pos}-{end_pos}'
            )
            chunk = response['Body'].read()

            newline_pos = chunk.index(b'\n')
            return (current_pos + newline_pos + 1, True)  # Found newline
        except ValueError:
            # No newline in this chunk, continue
            current_pos += chunk_size
            chunk_size = min(chunk_size * 2, 10*1024*1024)  # Double size, cap at 10MB
        except Exception as e:
            if 'InvalidRange' in str(e):
                break
            raise

    # No newline found in max_search range
    return (min(start_pos + max_search, file_size), False)  # Not found


def estimate_avg_line_length(s3, bucket, key, file_size, num_samples=5,
                             sample_size=256*1024, debug=False):
    """
    Estimate average line length by sampling multiple points across the file.

    Takes small samples from beginning, middle, end, and intermediate positions
    to get a representative average that handles variable line lengths.

    Args:
        s3: boto3 S3 client
        bucket: S3 bucket name
        key: S3 key
        file_size: Total file size in bytes
        num_samples: Number of sample points (default 5: start, 25%, 50%, 75%, end)
        sample_size: Bytes per sample (default 256KB)
        debug: Enable debug output

    Returns:
        Estimated average line length in bytes
    """
    try:
        total_bytes = 0
        total_newlines = 0

        # Calculate sample positions evenly distributed across file
        sample_positions = []
        if num_samples == 1:
            sample_positions = [0]
        else:
            for i in range(num_samples):
                # Distribute samples evenly: 0%, 25%, 50%, 75%, 100%
                position = int((file_size - sample_size) * i / (num_samples - 1))
                position = max(0, min(position, file_size - sample_size))
                sample_positions.append(position)

        if debug:
            print(f"[Sampling] File size: {file_size} bytes, taking {num_samples} samples of {sample_size} bytes each")

        for idx, start_pos in enumerate(sample_positions):
            end_pos = min(start_pos + sample_size - 1, file_size - 1)

            # Download sample
            response = s3.get_object(
                Bucket=bucket,
                Key=key,
                Range=f'bytes={start_pos}-{end_pos}'
            )
            sample = response['Body'].read()

            # Count newlines in this sample
            newline_count = sample.count(b'\n')
            total_newlines += newline_count
            total_bytes += len(sample)

            if debug:
                sample_avg = len(sample) / newline_count if newline_count > 0 else 0
                print(f"[Sampling] Sample {idx+1} @ {start_pos}: {newline_count} newlines, "
                      f"avg={sample_avg:.1f}B/line")

        if total_newlines == 0:
            # No newlines found - assume very long lines or binary data
            if debug:
                print(f"[Sampling] No newlines in {total_bytes} bytes across {num_samples} samples")
                print(f"[Sampling] Defaulting to 1MB window (likely binary or very long lines)")
            return 1 * 1024 * 1024  # Use 1MB window as safety

        avg_line_length = total_bytes / total_newlines

        if debug:
            print(f"[Sampling] Total: {total_bytes} bytes, {total_newlines} newlines")
            print(f"[Sampling] Estimated avg line length: {avg_line_length:.1f} bytes")

        return avg_line_length

    except Exception as e:
        if debug:
            print(f"[Sampling] Error during sampling: {e}")
            print(f"[Sampling] Defaulting to 64KB window")
        # Fallback to conservative default
        return 64 * 1024


def find_line_boundaries_smart(bucket, key, num_shards, chunks_per_lambda=1, window_size=None, debug=False):
    """
    Find line-aligned byte boundaries by downloading ONLY small windows around boundaries.

    If window_size is not provided, automatically determines optimal size by sampling file
    at multiple points (beginning, middle, end) to handle variable line lengths.

    Args:
        bucket: S3 bucket name
        key: S3 key
        num_shards: Number of lambda workers
        chunks_per_lambda: How many chunks each lambda processes
        window_size: Size of window to download around each boundary (None = auto-detect)
        debug: Enable debug output

    Returns:
        Tuple of (ranges, all_success) where:
        - ranges: List of (start_byte, end_byte, skip_first_line) tuples for each chunk
          - Has num_shards × chunks_per_lambda entries
          - skip_first_line=False if boundary was successfully line-aligned
          - skip_first_line=True if we had to fall back to approximate boundary
        - all_success: True if ALL boundaries were successfully found

    Example for 20GB file, 16 shards, 2 chunks per lambda (all successful):
        Total chunks: 16 × 2 = 32
        Downloads: 31 boundaries × 1MB = 31MB
        Returns: ([(0, 655359999, False), (655360000, 1310719999, False), ...], True)
    """
    with time_block("Total boundary scan"):
        s3 = boto3.client('s3')

        total_chunks = num_shards * chunks_per_lambda

        # Step 1: Get file size (needed for sampling positions)
        with time_block("S3 HEAD request"):
            response = s3.head_object(Bucket=bucket, Key=key)
            file_size = response['ContentLength']

        # Step 2: Adaptive window sizing
        if window_size is None:
            # Check for env var override first
            if 'PASH_BOUNDARY_WINDOW_KB' in os.environ:
                window_size_kb = int(os.environ.get('PASH_BOUNDARY_WINDOW_KB'))
                window_size = window_size_kb * 1024
                if debug:
                    print(f"[Boundary Scan] Using window_size={window_size/1024:.1f} KB (from env var)")
            else:
                # Multi-point sampling to estimate line length
                with time_block("Multi-point file sampling"):
                    # Take 5 samples of 256KB each = 1.25MB total
                    # Positions: 0%, 25%, 50%, 75%, 100% through file
                    avg_line_length = estimate_avg_line_length(
                        s3, bucket, key, file_size,
                        num_samples=5,
                        sample_size=256*1024,
                        debug=debug
                    )

                # Set window to hold ~500 lines (with 4KB minimum, 1MB maximum)
                # This provides safety margin while minimizing downloads
                target_lines_per_window = 500
                window_size = int(avg_line_length * target_lines_per_window)
                window_size = max(4 * 1024, min(window_size, 1024 * 1024))  # Clamp 4KB-1MB

                if debug:
                    print(f"[Boundary Scan] Adaptive window_size={window_size/1024:.1f} KB "
                          f"(avg_line={avg_line_length:.1f}B × {target_lines_per_window} lines)")

        all_boundaries_found = True  # Track if any boundary scan failed

        # Add verbosity level control
        verbose = debug and os.environ.get('PASH_VERBOSE_TIMING', 'false').lower() == 'true'

        if debug:
            print(f"[Boundary Scan] File size: {file_size} bytes, {num_shards} shards, {chunks_per_lambda} chunks/shard = {total_chunks} total chunks")

        # Edge case: Empty file
        if file_size == 0:
            if debug:
                print(f"[Boundary Scan] Empty file - returning trivial ranges")
            # All chunks get empty ranges, no skip needed
            return ([(0, 0, False) for _ in range(total_chunks)], True)

        # Edge case: Single chunk
        if total_chunks == 1:
            if debug:
                print(f"[Boundary Scan] Single chunk - no boundaries to scan")
            return ([(0, file_size - 1, False)], True)

        # Step 2: Calculate approximate chunk size
        approx_chunk_size = file_size // total_chunks

        # Track which boundaries were successfully found (vs fell back to approximate)
        # Index 0 is always 0 (start of file), so no skip needed for chunk 0
        boundary_success = [True]  # Chunk 0 starts at 0, always line-aligned

        boundaries = [0]  # First chunk starts at byte 0

        # Step 3: Find exact line boundaries
        with time_block(f"Boundary scan ({total_chunks-1} boundaries)"):
            for i in range(1, total_chunks):
                approx_boundary = i * approx_chunk_size
                found_newline = False

                # Download small window around the approximate boundary
                # 1MB window: 512KB before + 512KB after
                start = max(0, approx_boundary - window_size // 2)
                end = min(file_size - 1, approx_boundary + window_size // 2)

                if verbose:  # Only print per-chunk details if explicitly enabled
                    print(f"[Boundary Scan] Chunk {i}: downloading bytes {start}-{end} (window around {approx_boundary})")
                elif debug and i % 50 == 0:  # Print every 50th chunk as progress indicator
                    print(f"[Boundary Scan] Progress: {i}/{total_chunks} boundaries scanned")

                try:
                    response = s3.get_object(
                        Bucket=bucket,
                        Key=key,
                        Range=f'bytes={start}-{end}'
                    )
                    chunk = response['Body'].read()

                    # Find first newline AFTER the approximate boundary
                    offset_in_chunk = approx_boundary - start
                    try:
                        newline_pos = chunk.index(b'\n', offset_in_chunk)
                        # Boundary starts AFTER the newline (+1)
                        actual_boundary = start + newline_pos + 1
                        found_newline = True

                        if verbose:
                            print(f"[Boundary Scan] Chunk {i}: found newline at byte {actual_boundary}")

                        boundaries.append(actual_boundary)
                        boundary_success.append(True)

                    except ValueError:
                        # Newline not found in window - expand search
                        if verbose:
                            print(f"[Boundary Scan] Chunk {i}: newline not in window, expanding search...")

                        actual_boundary, found = expand_search_for_newline(
                            s3, bucket, key, approx_boundary, file_size
                        )
                        boundaries.append(actual_boundary)

                        if found:
                            boundary_success.append(True)
                            if verbose:
                                print(f"[Boundary Scan] Chunk {i}: found newline at byte {actual_boundary} (expanded search)")
                        else:
                            # No newline found even after expanding - use approximate, need skip
                            boundary_success.append(False)
                            all_boundaries_found = False
                            if verbose:
                                print(f"[Boundary Scan] Chunk {i}: NO NEWLINE FOUND, using approximate boundary {actual_boundary}, will skip first line")

                except Exception as e:
                    print(f"[Boundary Scan] ERROR at chunk {i}: {e}")
                    # Fallback to approximate boundary - will need skip_first_line
                    boundaries.append(approx_boundary)
                    boundary_success.append(False)
                    all_boundaries_found = False

        boundaries.append(file_size)  # Last chunk ends at EOF

        # Convert to (start, end, skip_first_line) ranges
        with time_block("Range calculation"):
            ranges = []
            for i in range(total_chunks):
                start_byte = boundaries[i]
                end_byte = boundaries[i + 1] - 1  # end_byte is inclusive

                # Validate range (handle edge case where boundaries collapse)
                if end_byte < start_byte:
                    # Empty or invalid range - give it at least one byte
                    end_byte = start_byte

                # skip_first_line is based on whether THIS chunk's START boundary was line-aligned
                # Chunk 0 starts at byte 0 (always line-aligned, no skip)
                # Chunk i>0 needs skip if boundary_success[i] is False
                skip_first_line = not boundary_success[i] if i > 0 else False

                ranges.append((start_byte, end_byte, skip_first_line))

        if debug:
            print(f"\n[Boundary Scan] Completed: {total_chunks} chunks, all_success={all_boundaries_found}")
            print(f"[Boundary Scan] Total data downloaded: ~{(total_chunks-1) * window_size / (1024*1024):.1f} MB")
            if verbose:  # Only print full range list if verbose mode enabled
                print(f"[Boundary Scan] Final ranges:")
                for i, (start, end, skip) in enumerate(ranges):
                    size_mb = (end - start + 1) / (1024 * 1024)
                    print(f"  Chunk {i}: bytes={start}-{end} ({size_mb:.2f} MB), skip={skip}")

        return (ranges, all_boundaries_found)

def get_s3_size(bucket, key):
    #print("Came here")
    """Get S3 object size."""
    s3 = boto3.client('s3')
    #print("bucket ", bucket)
    #print("key ", key)
    #print(boto3.client("sts").get_caller_identity())
    response = s3.head_object(Bucket=bucket, Key=key)
    #print("got obj")
    return response['ContentLength']

def add_nodes_to_subgraphs(subgraphs:List[IR], file_id_gen: FileIdGen, input_fifo_map:Dict[int, IR], args: argparse.Namespace, recover: bool = False):
    """ Takes a list of subgraphs and augments subgraphs with the necessary remote
        read/write nodes for data movement and lambda invocation nodes to trigger
        downstream processing (lambda). This function also produces graph that should run in
        the original shell in which pash was executed. EC2 (client graph) This graph contains
        remote read/write nodes for stdin/stdout, named pipes, and files.

    Args:
        subgraphs: list of sub sections of an optimized IR (returned from split_ir)
        file_id_gen: file id generator of the original ir
        input_fifo_map: mapping from input idge id to subgraph (returned from split_ir)
        args: command-line arguments including:
            - args.sls_output: S3 output location for final results
            - args.no_eager: Disable eager prefetching nodes
            - args.enable_s3_direct: Enable S3 direct lambda streaming optimization
        recover: whether to use recovery mode (default False)
    Returns:
        main_graph_script_id: the script id to execute on main shell
        subgraph_script_id_pairs: mapping from subgraph to unique script id
        main_subgraph_script_id: the script id to execute in the first lambda
    """
    # Print subgraphs for debugging/visualization

    if graph:
        print("\n" + "="*80)
        print("DEBUG: Visualizing subgraphs before adding nodes")
        print("="*80)

        pretty_print_subgraphs(subgraphs, unified_view=True)

        print("="*80)

    # The graph to execute in the main pash_compiler
    main_graph = IR({}, {})
    subgraph_script_id_pairs = {}
    main_subgraph_script_id = None

    subgraph_to_downstream_subgraphs = {} # subgraph -> [all downstream subgraphs], used for level order traversal
    subgraph_to_invocation_node = {} # subgraph -> the invocation node invoking this subgraph
    subgraphs_not_having_upstream = set(subgraphs)
    subgraph_stun_lib_args = {}
    key_to_sender_receiver = {}
    key_to_data_type = {} # batch or line
    eager_edges = []

    #modify the first subgraph if it is a s3 -> cat -> split
    # then the loop will not even have it anymore
    # helper fn
    #subgraphs = change_sorting(subgraphs)



    ec2_in_edge = None
    total_lambdas = 0
    if args.enable_s3_direct:
        print("[IR Helper] S3 direct streaming optimization ENABLED")
        subgraphs, ec2_in_edge, total_lambdas = optimize_s3_lambda_direct_streaming(subgraphs, input_fifo_map)
    else:
        print("[IR Helper] S3 direct streaming optimization DISABLED (use --enable_s3_direct to enable)")
        # Keep original subgraphs unchanged - matching the "no optimization" return pattern

    #generate graphviz here so we can see the diff
    # actually not trivial as we have a list of IRs here :( not a single IR 

    print("AFTER OPTIMIZATION, TOTAL LAMBDAS STREAMING DIRECTLY FROM S3:", total_lambdas)

    if graph:
        print("\n" + "="*80)
        print("DEBUG: Visualizing subgraphs after S3 optimization")
        print("="*80)

        pretty_print_subgraphs(subgraphs, unified_view=True)

        print("="*80)
    #exit()
    #TODO careful with this exit
    
    # Replace output edges and corrosponding input edges with remote read/write
    # with the key as old_edge_id
    for subgraph in subgraphs:
        sink_nodes = subgraph.sink_nodes()
        assert(len(sink_nodes) == 1)
        out_edges = subgraph.get_node_output_fids(sink_nodes[0])
        for out_edge in out_edges:
            # Replace the old edge with an ephemeral edge in case it isn't and
            # to avoid modifying the edge in case it's used in some other subgraph
            out_edge_id = out_edge.get_ident()
            ephemeral_edge = file_id_gen.next_ephemeral_file_id()
            subgraph.replace_edge(out_edge_id, ephemeral_edge)
            communication_key = uuid4()
            stdout = add_stdout_fid(subgraph, file_id_gen)
            # if no downstream subgraph, assuming this is the last subgraph
            # TODO: check if above assumption makes sense
            last_subgraph = False
            if out_edge_id not in input_fifo_map:
                last_subgraph = True
            if last_subgraph:
                communication_key = "stdout"
                if type(out_edge.get_resource()) is FileResource:
                    communication_key = str(out_edge.get_resource())
                if args.sls_output != "":
                    communication_key = os.path.join(args.sls_output, str(communication_key))
            # Add remote-write node at the end of the subgraph
            if (not last_subgraph):
                # arg = send_rdvkey_0_1_input-fifoname
                out_node = subgraph.get_node(sink_nodes[0])
                # print(out_node)
                key_to_data_type[str(communication_key)] = "line"
                if isinstance(out_node, RWrap):
                    key_to_data_type[str(communication_key)] = "batch"
                if isinstance(out_node, RSplit):
                    key_to_data_type[str(communication_key)] = "batch"
                    for flag in out_node.cmd_invocation_with_io_vars.flag_option_list:
                        if flag.get_name() == "-r":
                            key_to_data_type[str(communication_key)] = "line"
                if str(communication_key) not in key_to_sender_receiver:
                    key_to_sender_receiver[str(communication_key)] = [subgraph, None]
                else:
                    key_to_sender_receiver[str(communication_key)][1] = subgraph
                
                # sends lambda's output to ec2/lambda 
                # not sure. have to figure that out
                arg = "send*"+str(communication_key)+"*0*1*"+config.PASH_TMP_PREFIX+str(ephemeral_edge) 
                
                if subgraph not in subgraph_stun_lib_args:
                    subgraph_stun_lib_args[subgraph] = []
                subgraph_stun_lib_args[subgraph].append(arg)
            else:
                remote_write = serverless_remote_pipe.make_serverless_remote_pipe(local_fifo_id=ephemeral_edge.get_ident(),
                                                                                is_remote_read=False,
                                                                                remote_key=communication_key,
                                                                                output_edge=stdout,
                                                                                is_tcp=(not last_subgraph))
                subgraph.add_node(remote_write)

            # Copy the old output edge resource
            new_edge = file_id_gen.next_file_id()
            new_edge.set_resource(out_edge.get_resource())
            # Get the subgraph which "edge" writes to
            if out_edge_id in input_fifo_map and out_edge.is_ephemeral():
                matching_subgraph = input_fifo_map[out_edge_id][0]
                matching_subgraph.replace_edge(out_edge.get_ident(), new_edge)
                # Add invocation node
                if matching_subgraph not in subgraph_script_id_pairs:
                    script_identifier = uuid4()
                    subgraph_script_id_pairs[matching_subgraph] = script_identifier
                    # incovation_node = serverless_lambda_invoke.make_serverless_lambda_invoke(script_identifier)
                    # subgraph.add_node(incovation_node)
                    # subgraph_to_invocation_node[matching_subgraph] = incovation_node
                if subgraph not in subgraph_to_downstream_subgraphs:
                    subgraph_to_downstream_subgraphs[subgraph] = []
                # print("Adding downstream subgraph to ", subgraph_script_id_pairs[matching_subgraph])
                subgraph_to_downstream_subgraphs[subgraph].append(matching_subgraph)
                subgraphs_not_having_upstream.discard(matching_subgraph)
            else:
                # Add edge to main graph
                #TODO. inspect subgraph obj to see before and after 
                # defined in compiler/IR.py

                matching_subgraph = main_graph
                matching_subgraph.add_edge(new_edge)
            if (not matching_subgraph is main_graph):
                # arg = recv_rdvkey_1_0_output-fifoname

                # remote rcv in ec2/lambda from ec2/lambda
                arg = "recv*"+str(communication_key)+"*1*0*"+config.PASH_TMP_PREFIX+str(new_edge)
                if matching_subgraph not in subgraph_stun_lib_args:
                    subgraph_stun_lib_args[matching_subgraph] = []
                subgraph_stun_lib_args[matching_subgraph].append(arg)

                if str(communication_key) not in key_to_sender_receiver:
                    key_to_sender_receiver[str(communication_key)] = [None, matching_subgraph]
                else:
                    key_to_sender_receiver[str(communication_key)][1] = matching_subgraph

                if not args.no_eager:
                    if recover:
                        eager_edges.append((new_edge, matching_subgraph))
                    else:
                        pash_compiler.add_eager(new_edge.get_ident(), matching_subgraph, file_id_gen)

            else: #similar to what is done here we want to add this node before the lambda pash node or in lieu of
                remote_read = serverless_remote_pipe.make_serverless_remote_pipe(local_fifo_id=new_edge.get_ident(),
                                                                                is_remote_read=True,
                                                                                remote_key=communication_key,
                                                                                output_edge=new_edge,
                                                                                is_tcp=(not matching_subgraph is main_graph))
                matching_subgraph.add_node(remote_read)

    if graph:
        print("\n\nAFTER FIRST FOR LOOP")
        print("="*80)

        pretty_print_subgraphs(subgraphs, unified_view=True)

        print("="*80)

    lambda_counter = 0
    job_uid = uuid4()

    shard_ranges = None  # Will be populated on first S3 lambda
    shard_ranges_computed = False
    # Replace non ephemeral input edges with remote read/write
    for subgraph in subgraphs:
        if subgraph not in subgraph_script_id_pairs:
            main_subgraph_script_id = uuid4()
            subgraph_script_id_pairs[subgraph] = main_subgraph_script_id
        source_nodes = subgraph.source_nodes()
        for source in source_nodes:
            if isinstance(subgraph.get_node(source), serverless_lambda_invoke.ServerlessLambdaInvoke) or \
                isinstance(subgraph.get_node(source), serverless_remote_pipe.ServerlessRemotePipe):
                continue
            for in_edge in subgraph.get_node_input_fids(source):
                # TODO: also consider in_edge.has_file_descriptor_resource()
                if in_edge.has_file_resource():
                    filename = in_edge.get_resource().uri
                    # Add remote read to current subgraph
                    ephemeral_edge = file_id_gen.next_ephemeral_file_id()
                    subgraph.replace_edge(in_edge.get_ident(), ephemeral_edge)
                    

                    if in_edge == ec2_in_edge:  # TODO how to get bucket??
                        BUCKET=os.environ.get("AWS_BUCKET")
                        filesize = get_s3_size(BUCKET, str(filename).strip('"'))
                        # but filename can be covid-mts/inputs/in.csv so we need a more robust way
                        # Extract file size from filename (e.g., "oneliners/inputs/1M.txt" -> 1048576)
                        # size_multipliers = {"G": 1024**3, "M": 1024**2, "K": 1024}
                        filename_stripped = str(filename).strip('"')
                        print(f"File: {filename} -> Size: {filesize} bytes")

                        # Read configuration
                        chunks_per_lambda = int(os.environ.get('PASH_S3_CHUNKS_PER_LAMBDA', '1'))
                        if chunks_per_lambda < 1:
                            chunks_per_lambda = 1
                            print(f"[IR Helper] WARNING: Invalid PASH_S3_CHUNKS_PER_LAMBDA, using 1")

                        if chunks_per_lambda > 1:
                            print(f"[IR Helper] Multi-chunk mode: {chunks_per_lambda} chunks per lambda")
                        else:
                            print(f"[IR Helper] Single-chunk mode (original behavior)")

                        # Check if smart boundaries are enabled
                        use_smart_boundaries = os.environ.get('USE_SMART_BOUNDARIES', 'false').lower() == 'true'

                        if use_smart_boundaries:
                            # Calculate line-aligned boundaries (only once for all lambdas)
                            # Note: This will be called once per lambda, but we only need to calculate once
                            # For now, we calculate per lambda (could be optimized to cache)
                          

                            if not shard_ranges_computed:
                                print(f"[IR Helper] Scanning S3 file for line boundaries: {filename_stripped}")
                                print(f"[IR Helper] Total lambdas: {total_lambdas}, chunks per lambda: {chunks_per_lambda}")

                                with time_block("find_line_boundaries_smart call"):
                                    shard_ranges, boundary_scan_success = find_line_boundaries_smart(
                                        bucket=BUCKET,
                                        key=filename_stripped,
                                        num_shards=total_lambdas,
                                        chunks_per_lambda=chunks_per_lambda,
                                        window_size=None,  # Enable adaptive sizing
                                        debug=True
                                    )
                                shard_ranges_computed = True

                                if not boundary_scan_success:
                                    print(f"[IR Helper] WARNING: Boundary scan had failures")

                            # Assign chunks to THIS lambda (round-robin interleaved)
                            with time_block(f"Chunk distribution (lambda {lambda_counter})"):
                                lambda_chunks = []
                                for chunk_index in range(chunks_per_lambda):
                                    global_chunk_id = lambda_counter + (chunk_index * total_lambdas)

                                    # Safety check
                                    if global_chunk_id >= len(shard_ranges):
                                        break

                                    start_byte, end_byte, use_skip = shard_ranges[global_chunk_id]
                                    lambda_chunks.append({
                                        "start": start_byte,
                                        "end": end_byte,
                                        "block_id": global_chunk_id,  # CRITICAL: Global sequential ID
                                        "skip_first_line": use_skip
                                    })

                            # Format byte_range parameter based on mode
                            if chunks_per_lambda == 1:
                                # Single-chunk mode: backward compatible string format
                                chunk = lambda_chunks[0]
                                byte_range = f"bytes={chunk['start']}-{chunk['end']}"
                                skip_first_line = chunk['skip_first_line']

                                print(f"[IR Helper] Lambda {lambda_counter}: {byte_range}, skip_first_line={skip_first_line}")
                            else:
                                # Multi-chunk mode: JSON array
                                import shlex
                                byte_range_json = json.dumps(lambda_chunks)                                                             
                                byte_range = f"'{byte_range_json}'"
                                # byte_range = json.dumps(lambda_chunks)
                                skip_first_line = False  # Not used in multi-chunk mode

                                total_mb = sum((c['end'] - c['start'] + 1) / (1024*1024) for c in lambda_chunks)
                                chunk_ids = [c['block_id'] for c in lambda_chunks]
                                print(f"[IR Helper] Lambda {lambda_counter}: {len(lambda_chunks)} chunks, IDs {chunk_ids}, {total_mb:.2f} MB total")


                        else:
                            # Use approximate boundaries (old behavior)
                            skip_first_line = True
                            chunk_size = filesize // total_lambdas
                            start_byte = lambda_counter * chunk_size
                            # Last lambda gets everything remaining to handle rounding
                            end_byte = filesize - 1 if lambda_counter == total_lambdas - 1 else (lambda_counter + 1) * chunk_size - 1

                            byte_range = f"bytes={start_byte}-{end_byte}"
                            print(f"Lambda {lambda_counter}: {byte_range}")
                        remote_read = serverless_remote_pipe.make_serverless_remote_pipe(local_fifo_id=ephemeral_edge.get_ident(),
                                                                                is_remote_read=True,
                                                                                remote_key=filename,
                                                                                output_edge=None,
                                                                                is_tcp=False,
                                                                                is_s3_lambda=True,
                                                                                lambda_counter=lambda_counter,
                                                                                total_lambdas=total_lambdas,
                                                                                byte_range=byte_range, 
                                                                                job_uid=job_uid,
                                                                                skip_first_line=skip_first_line)
                        
                        # s3getobj byterange=start-end -> sort   : lambda1
                        lambda_counter += 1
                    else:
                        remote_read = serverless_remote_pipe.make_serverless_remote_pipe(local_fifo_id=ephemeral_edge.get_ident(),
                                                                                is_remote_read=True,
                                                                                remote_key=filename,
                                                                                output_edge=None,
                                                                                is_tcp=False)
                    if in_edge == ec2_in_edge and not args.no_eager:
                        # Add dgsh-tee for eager S3 data prefetching when using S3 direct streaming
                        # This ensures data is pulled from S3 as fast as possible and buffered for downstream
                
                        pash_compiler.add_eager(ephemeral_edge.get_ident(), subgraph, file_id_gen, is_s3=True)

                    # print("Adding remote read for S3 file:", filename)
                    subgraph.add_node(remote_read) # This makes a node with an s3 input (converts the input filename to an s3 get obj cmd)
                    # print("Remote read node added:", remote_read)
                    

                else:
                    # sometimes a command can have both a file resource and an ephemeral resources (example: spell oneliner)
                    continue

    if graph:
        print("\n\nAFTER SECOND FOR LOOP (where we add s3 read nodes, i.e convert the input edges into s3 get nodes)")
        print("="*80)

        pretty_print_subgraphs(subgraphs, unified_view=True)

        print("="*80)

    # Print timing summary if boundary scanning was performed
    if shard_ranges_computed:
        print_timing_summary()

    main_graph_script_id = uuid4()
    subgraph_script_id_pairs[main_graph] = main_graph_script_id

    fifo_to_be_renamed = {}
    if recover:
        subgraph_types = {}
        lambda_upstream_downstream = {} #?
        for subgraph, args in subgraph_stun_lib_args.items():
            send_count = recv_count = 0
            for arg in args:
                if "send" in arg:
                    send_count += 1
                else:
                    recv_count += 1
            if send_count == 1 and recv_count == 1:
                subgraph_types[subgraph] = "lambda"
                sender_key = ""
                recver_key = ""
                for arg in args:
                    if "send" in arg:
                        sender_key = arg.split("*")[1]
                    else:
                        recver_key = arg.split("*")[1]
                downstream_subgraph = key_to_sender_receiver[sender_key][1]
                upstream_subgraph = key_to_sender_receiver[recver_key][0]
                lambda_upstream_downstream[subgraph] = (upstream_subgraph, downstream_subgraph)
            elif send_count > 1 and recv_count > 1:
                subgraph_types[subgraph] = "splitter_merger"
            elif send_count > 1:
                subgraph_types[subgraph] = "splitter"
            elif recv_count > 1:
                subgraph_types[subgraph] = "merger"
            else:
                subgraph_types[subgraph] = "unknown"

        # add eager for all non-lambda
        for edge, subgraph in eager_edges:
            if subgraph_types[subgraph] != "lambda":
                pash_compiler.add_eager(edge.get_ident(), subgraph, file_id_gen)
        
        # add ingate outgate
        for subgraph, stun_lib_args in subgraph_stun_lib_args.items():
            if subgraph_types[subgraph] == "lambda":
                # change all key to key_v${version}
                for i in range(len(stun_lib_args)):
                    key = stun_lib_args[i].split("*")[1]
                    stun_lib_args[i] = stun_lib_args[i].replace(key, key+"_v0")
            else:
                # for all send, change key to key_v0, change fifo to fifo_v0
                for i in range(len(stun_lib_args)):
                    key = stun_lib_args[i].split("*")[1]
                    stun_lib_args[i] = stun_lib_args[i].replace(key, key+"_v0")
                    fifo = stun_lib_args[i].split("*")[4]
                    stun_lib_args[i] = stun_lib_args[i].replace(fifo, fifo+"_v0")
                    # fifo_to_be_renamed.add(fifo)
                    if "send" in stun_lib_args[i]:
                        receiver = key_to_sender_receiver[key][1]
                        receiver_script_id = subgraph_script_id_pairs[receiver]
                        outgate_args = [fifo, fifo, key, "900", str(receiver_script_id)]
                        outgate = serverless_remote_pipe.make_serverless_outgate(key_to_data_type[key], outgate_args)
                        subgraph.add_node(outgate)
                        new_fid = file_id_gen.next_ephemeral_file_id()
                        subgraph.add_edge(new_fid)
                        if subgraph not in fifo_to_be_renamed:
                            fifo_to_be_renamed[subgraph] = set()
                        fifo_to_be_renamed[subgraph].add((str(config.PASH_TMP_PREFIX)+str((new_fid.get_fifo_suffix())), fifo+"_v0"))
                    else:
                        sender = key_to_sender_receiver[key][0]
                        sender_script_id = subgraph_script_id_pairs[sender]
                        ingate_args = [fifo, key, fifo, str(sender_script_id)]
                        ingate = serverless_remote_pipe.make_serverless_ingate(key_to_data_type[key], ingate_args)
                        subgraph.add_node(ingate)
                        new_fid = file_id_gen.next_ephemeral_file_id()
                        subgraph.add_edge(new_fid)
                        if subgraph not in fifo_to_be_renamed:
                            fifo_to_be_renamed[subgraph] = set()
                        fifo_to_be_renamed[subgraph].add((str(config.PASH_TMP_PREFIX)+str((new_fid.get_fifo_suffix())), fifo+"_v0"))
            subgraph_stun_lib_args[subgraph] = stun_lib_args

    for subgraph, stun_lib_args in subgraph_stun_lib_args.items():
        # stun_lib = serverless_remote_pipe.make_serverless_remote_pipe_one_proc(stun_lib_args)
        # subgraph.add_node(stun_lib)

        args_lists = [[]]
        peer_lists = [set()]
        for arg in stun_lib_args:
            if recover:
                key = arg.split("*")[1].split("_v")[0]
            else:
                key = arg.split("*")[1]
            sender, receiver = key_to_sender_receiver[key]
            if "send" in arg:
                peer = receiver
            else:
                peer = sender
            added = False
            for i in range(len(args_lists)):
                if peer not in peer_lists[i]:
                    args_lists[i].append(arg)
                    peer_lists[i].add(peer)
                    added = True
                    break
            if not added:
                args_lists.append([arg])
                peer_lists.append({peer})
        for args_list in args_lists:
            stun_lib = serverless_remote_pipe.make_serverless_remote_pipe_one_proc(args_list)
            subgraph.add_node(stun_lib)

    return main_graph_script_id, subgraph_script_id_pairs, main_subgraph_script_id, fifo_to_be_renamed


def prepare_scripts_for_serverless_exec(ir: IR, shell_vars: dict, args: argparse.Namespace, declared_functions_filename, recover: bool = False) -> Tuple[str, str, Dict[str, str]]:
    """
    Reads the complete ir from filename and splits it
    into subgraphs where ony the first subgraph represent a continues
    segment (merger segment or branched segment) in the graph.
    Note: All subgraphs(except first one) read and write from remote pipes.
        However, we had to add a fake stdout to avoid some problems when converting to shell code.

    Args:
        ir: an IR optimized by the pash compiler
        shell_vars: shell variables
        args: pash args

    Returns:
        main_graph_script_id: the script id to execute on main shell
        main_subgraph_script_id: the script id to execute in the first lambda
        script_id_to_script: mapping from unique script id to script content
    """
    # split IR
    subgraphs, mapping = split_ir(ir)

    #todo. chaneg the first one 
    main_graph_script_id, subgraph_script_id_pairs, main_subgraph_script_id, fifo_to_be_replaced = add_nodes_to_subgraphs(subgraphs, ir.get_file_id_gen(), mapping, args, recover=recover)

    # read the declared functions
    declared_functions = ""
    with open(declared_functions_filename, "r") as f:
        declared_functions = f.read()
    
    # save the output scripts
    script_id_to_script = {}
    ec2_set = set()
    for subgraph, id_ in subgraph_script_id_pairs.items():
        # making necessary temp directories
        dir_set = set()
        for edge in subgraph.all_fids():
            if edge.is_ephemeral():
                dir_set.add(os.path.join(config.PASH_TMP_PREFIX, edge.prefix))
        mk_dirs = "mkdir -p "+config.PASH_TMP_PREFIX+" \n"
        for dir in dir_set:
            mk_dirs += "mkdir -p "+dir+" \n"
        # rust_debug = "export RUST_BACKTRACE=full\n"
        export_path = "export PATH=$PATH:runtime\n"
        export_rust_trace = "export RUST_BACKTRACE=1\n"
        export_lib_path = "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:runtime/lib\n"
        export_locale_path = "export LOCPATH=/var/task/runtime/locale\n" 
        export_lang = "export LANG=C.UTF-8\n"
        export_locale_all = "export LOCALE_ALL=C.UTF-8\n"
        add_version = "version=$2\n"
        script = export_path+export_lib_path+export_locale_path+export_lang+export_locale_all+export_rust_trace+add_version+mk_dirs+f"{declared_functions}\n"+to_shell(subgraph, args)
        # generate scripts
        if recover and subgraph in fifo_to_be_replaced:
            for new_fifo, recover_fifo in fifo_to_be_replaced[subgraph]:
                script = script.replace(new_fifo, recover_fifo)
        script_name = os.path.join(config.PASH_TMP_PREFIX, str(id_))
        script_id_to_script[str(id_)] = script
        with open (script_name, "w") as f:
            f.write(script)
        if id_ == main_graph_script_id:
            log("[Serverless Manager] Script for main shell saved in:"+script_name)
        elif id_ == main_subgraph_script_id:
            log("[Serverless Manager] Script for first lambda saved in:"+script_name)
        else:
            log("[Serverless Manager] Script for other lambda saved in:"+script_name)
        # log(script)
        if ("split" in script) or ("s3-put" in script) or ("sort -m" in script) or ("merge" in script):
            ec2_set.add(str(id_))

    return str(main_graph_script_id), str(main_subgraph_script_id), script_id_to_script, ec2_set
