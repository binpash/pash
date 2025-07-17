# dspash – FRACTAL Compiler / Scheduler Extensions

**FRACTAL** brings *fault-tolerant* distribution to PaSh-JIT.  The code in this
directory implements the control-plane described in our (conditionally accepted)
NSDI’26 paper “Fractal: Fault-Tolerant Shell-Script Distribution” (citation
forthcoming).

```
compiler/dspash/
├── ir_helper.py        # subgraph carving, EdgeID generation, remote-pipe insertion
├── worker_manager.py   # global scheduler, batching, fault handling
├── worker.py           # per-node executor event-loop (B1-4)
├── hdfs_utils.py       # JMX polling for Health Monitor (A6)
└── __init__.py
```

## Build-time vs. Run-time
PaSh-JIT still performs JIT expansion & parallelization.  **dspash** kicks in
*after* a data-flow graph (DFG) is generated:
1. `prepare_graph_for_remote_exec` isolates the *unsafe main* region and tags
   the remaining vertices with HDFS block information.
2. `split_main_graph` divides the graph into *regular*, *merger*, and *singular*
   subgraphs and assigns a globally unique `SubgraphID`.
3. `add_remote_pipes` inserts proxy vertices and allocates persistent buffers
   when dynamic-persistence is on.  (Implementation is split across
   `create_remote_pipe_from_output_edge` / `create_remote_pipe_from_input_edge`
   in `ir_helper.py`.)

The resulting per-subgraph shell snippets are serialised (JSON) and shipped to
remote **worker** processes.

## Scheduler & Batching (worker_manager.py)
The coordinator tracks state in dictionaries (`all_worker_subgraph_pairs`,
`all_graph_to_uuid`, etc. – see `worker_manager.py`).

Its main loop reacts to: (1) new execution requests, (2) 17-byte completion
events from executors, and (3) health-monitor callbacks.  When batching is
enabled it groups subgraphs that share the same destination worker into a
single gRPC call.

Executor event-loops poll every **100 ms** (`worker.py`, `EventLoop.quit.wait(0.1)`).

Specifically the manager:
1. **Batches** subgraphs targeting the same node into a single gRPC request.
2. Issues `KillSubgraph` RPCs for anything that became unreachable per Health
   Monitor.
3. Checks `completion_events` to trigger *minimal* re-execution after faults.

## Executor Event Loop (`worker.py`)
A lock-free loop (*poll interval = 100 ms*) that:
1. Launches up to `concurrency = $CPU` sub-processes.
2. Reaps finished children and appends 17-byte events to the Progress Monitor.
3. Garbage-collects persisted files when both ends of an edge completed.

This replaces heavyweight container orchestration found in Hadoop.

## Health Monitor (hdfs_utils.py)
Queries the Namenode’s `/jmx` endpoint for `LastContact` and flags nodes whose
heartbeat > 10s, consistent to the default in HDFS.

## Dynamic Output Persistence
FRACTAL can decide **at runtime, per-edge**, whether writing to disk is cheaper
than replaying computation on a failure.  The decision path is:

1. When `--ft dynamic` is used, the scheduler tags *every* subgraph
   with the system-wide FT mode `dynamic` (see the `--ft` flag injected in
   `make_remote_pipe()` in `remote_pipe.py`).
2. `ir_helper.add_singular_flags()` labels subgraphs that have no downstream
   dependants ("Fig. 3 singular subgraphs"); their `RemotePipe` commands gain
   the `-s` flag so they **never** persist.
3. During execution the `datastream` client evaluates

   ```go
   useOptimizedRW := *ft == "optimized" || (*ft == "dynamic" && !*singular)
   ```

   • `optimized` → always persist  
   • `dynamic`   → persist unless `-s`.

4. After a fault the coordinator calls `worker_manager.check_persisted_discovery()`
   to ask Discovery which streams already exist on stable storage and thereby
   **skips** re-executing upstream subgraphs whose outputs have been cached.

This fine-grained path embodies the “dynamic output persistence” knob discussed
in §4.2 / Fig. 3-A3 of the paper.

## Relation to Paper
Implements Fig. 3-A1 (DFG augmentation), A3 (dynamic persistence policy),
Fig. 3-A4 (scheduler & batching), Fig. 3-A5 (progress monitor plumbing) and
Fig. 3-A6 (health monitor). 