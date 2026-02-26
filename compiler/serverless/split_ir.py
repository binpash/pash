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
from serverless.graph_print_helper import pretty_print_subgraphs
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
from definitions.ir.dfg_node import DFGNode
from definitions.ir.arg import Arg
from collections import defaultdict, deque
from typing import Dict, List, Tuple
from collections import defaultdict, deque
from typing import Dict, List, Tuple

def split_ir(graph: IR) -> Tuple[List[IR], Dict[int, IR]]:
    """
    Split only at the OUTERMOST split-merge pair:
      - split at the FIRST fork node (split_id)
      - merge at the FIRST reconvergence node of those branches (merge_id)
    Returns: [prefix] + [one per branch] + [suffix]  (=> N + 2 when merge exists)
    """

    # ---------------- helpers to find outermost split/merge ----------------
    def _collect_reachable_nodes(g: IR) -> List[int]:
        seen = set()
        q = deque(g.source_nodes())
        while q:
            u = q.popleft()
            if u in seen:
                continue
            seen.add(u)
            for v in g.get_next_nodes(u):
                if v not in seen:
                    q.append(v)
        return list(seen)

    def _topo_order_reachable(g: IR) -> List[int]:
        nodes = _collect_reachable_nodes(g)
        node_set = set(nodes)
        indeg = {nid: 0 for nid in nodes}
        for u in nodes:
            for v in g.get_next_nodes(u):
                if v in node_set:
                    indeg[v] += 1

        q = deque([nid for nid in nodes if indeg[nid] == 0])
        order = []
        while q:
            u = q.popleft()
            order.append(u)
            for v in g.get_next_nodes(u):
                if v in indeg:
                    indeg[v] -= 1
                    if indeg[v] == 0:
                        q.append(v)
        return order if len(order) == len(nodes) else nodes

    def _reachable_from(g: IR, start_id: int) -> set:
        seen = set()
        stack = [start_id]
        while stack:
            u = stack.pop()
            if u in seen:
                continue
            seen.add(u)
            stack.extend(g.get_next_nodes(u))
        return seen

    def _find_outermost_split_merge(g: IR):
        topo = _topo_order_reachable(g)
        topo_idx = {nid: i for i, nid in enumerate(topo)}

        # first fork in topo order
        split_id = None
        for nid in topo:
            if len(g.get_next_nodes(nid)) > 1:
                split_id = nid
                break
        if split_id is None:
            return None, None

        succs = g.get_next_nodes(split_id)
        if not succs:
            return split_id, None

        common = None
        for s in succs:
            r = _reachable_from(g, s)
            common = r if common is None else (common & r)
        if not common:
            return split_id, None

        def indegree(nid: int) -> int:
            return len(g.get_node_input_fids(nid))

        # earliest common node with indegree > 1
        candidates = [nid for nid in common if indegree(nid) > 1]
        if not candidates:
            return split_id, None

        merge_id = min(candidates, key=lambda nid: topo_idx.get(nid, 10**18))
        return split_id, merge_id

    split_id, merge_id = _find_outermost_split_merge(graph)

    # ---------------- splitting logic ----------------
    source_node_ids = graph.source_nodes()
    input_fifo_map = defaultdict(list)

    subgraphs: List[IR] = []
    queue = deque([(source, IR({}, {})) for source in source_node_ids])

    # keep track of traversed edges (readiness)
    visited_edges = set(graph.all_input_fids())

    # prevent re-processing the same (node, subgraph_object) pair
    visited = set()

    # prevent appending the SAME IR object multiple times
    finalized = set()

    def flush(sg: IR):
        if sg.source_nodes() and id(sg) not in finalized:
            subgraphs.append(sg)
            finalized.add(id(sg))

    did_outer_split = False
    started_suffix = False

    while queue:
        old_node_id, subgraph = queue.popleft()

        # If we're inside a branch subgraph and we reached the OUTER merge,
        # close the branch *before* consuming merge, and start suffix once.
        if (
            merge_id is not None
            and did_outer_split
            and old_node_id == merge_id
            and subgraph.source_nodes()
        ):
            flush(subgraph)  # keep distinct branch IRs; avoid duplicate same-object appends
            if not started_suffix:
                started_suffix = True
                queue.appendleft((old_node_id, IR({}, {})))  # re-process merge as start of suffix
            continue

        key = (old_node_id, id(subgraph))
        if key in visited:
            continue
        visited.add(key)

        input_fids = graph.get_node_input_fids(old_node_id)
        output_fids = graph.get_node_output_fids(old_node_id)

        # Not ready: end current linear growth (if any), but DO NOT double-append.
        if any(fid not in visited_edges for fid in input_fids):
            flush(subgraph)
            continue

        node = graph.get_node(old_node_id).copy()
        node_id = node.get_id()

        # Attach incoming edges and record input_fifo_map
        for input_fid in input_fids:
            if input_fid.get_ident() not in subgraph.edges:
                subgraph.add_to_edge(input_fid, node_id)
                input_edge_id = input_fid.get_ident()
            else:
                input_edge_id = input_fid.get_ident()
                subgraph.set_edge_to(input_edge_id, node_id)
            input_fifo_map[input_edge_id].append(subgraph)

        # Add outgoing edges
        for output_fid in output_fids:
            subgraph.add_from_edge(node_id, output_fid)
            visited_edges.add(output_fid)

        # Add the node
        subgraph.add_node(node)

        next_ids = graph.get_next_nodes(old_node_id)

        # Sink
        if len(next_ids) == 0:
            flush(subgraph)
            continue

        # Linear
        if len(next_ids) == 1:
            queue.append((next_ids[0], subgraph))
            continue

        # Fork: ONLY cut at the OUTER split_id (the "matching split")
        if split_id is not None and old_node_id == split_id and not did_outer_split:
            did_outer_split = True
            # close prefix (includes the split node)
            flush(subgraph)
            # start each branch as fresh subgraph
            for next_id in next_ids:
                queue.append((next_id, IR({}, {})))
        else:
            # nested forks: do NOT split; keep going with same subgraph
            for next_id in next_ids:
                queue.append((next_id, subgraph))

    return subgraphs, input_fifo_map