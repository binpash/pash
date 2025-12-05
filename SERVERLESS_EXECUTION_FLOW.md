# PaSh Serverless Execution Flow - Complete Function Call Chain

**Created:** 2025-11-25
**Purpose:** Deep dive into how ServerlessManager, ir_helper, and graph construction work together

---

## Overview: The Complete Flow

```
User Request
    ↓
ServerlessManager.handler()
    ↓
prepare_scripts_for_serverless_exec()
    ↓
split_ir() → add_nodes_to_subgraphs() → to_shell()
    ↓
Upload to S3 & Invoke Lambda/EC2
```

---

## 1. Entry Point: ServerlessManager.handler()

**File:** `/home/ubuntu/pash/compiler/serverless/serverless_executor.py`

**Location:** Lines 213-294

### What It Does

The `handler()` method is called when a request arrives via socket. It:
1. Parses the request to get IR filename and declared functions file
2. Calls `init()` to load the IR graph
3. Calls `prepare_scripts_for_serverless_exec()` to generate scripts
4. Uploads scripts to S3
5. Invokes Lambda/EC2(local) workers

### Code Flow

```python
def handler(self, request, conn):
    # 1. Parse request
    args = request.split(':', 1)[1].strip()
    ir_filename, declared_functions_file = args.split()

    # 2. Load IR graph from pickled file
    ir, shell_vars, args = init(ir_filename)

    # 3. MAIN CALL: Generate scripts for serverless execution
    main_graph_script_id, main_subgraph_script_id, script_id_to_script, ec2_set = \
        prepare_scripts_for_serverless_exec(
            ir,                          # The optimized IR graph
            shell_vars,                  # Shell variables
            args,                        # PaSh arguments
            declared_functions_file,     # Bash functions to include
            recover=self.recovery        # Recovery mode flag
        )

    # 4. Generate unique folder ID for this execution
    s3_folder_id = str(int(time.time()))  # Timestamp

    # 5. Upload all scripts to S3 (multi-threaded)
    s3 = boto3.client('s3', config=client_config)
    bucket = os.getenv("AWS_BUCKET")

    for script_id, script in script_id_to_script.items():
        if script_id == main_graph_script_id:
            continue  # Skip main script
        # Upload: s3://bucket/sls-scripts/{timestamp}/{script_id}.sh
        s3.put_object(
            Bucket=bucket,
            Key=f'sls-scripts/{s3_folder_id}/{script_id}.sh',
            Body=script
        )

    # 6. Respond to caller with folder ID
    response_msg = f"OK {s3_folder_id} {main_subgraph_script_id}"
    conn.sendall(response_msg.encode('utf-8'))

    # 7. Invoke all workers (EC2 or Lambda)
    for script_id, script in script_id_to_script.items():
        if script_id == main_graph_script_id:
            continue

        if self.ec2_enabled and script_id in ec2_set:
            # Run on EC2 (locally)
            invocation_thread = threading.Thread(
                target=self.run_local,
                args=(script_id,)
            )
        else:
            # Run on Lambda
            invocation_thread = threading.Thread(
                target=self.invoke_lambda,
                args=([s3_folder_id], [script_id])
            )
        invocation_thread.start()
```

### Key Insight

The `handler()` method is the orchestrator. It doesn't understand the details of how scripts are generated - it just calls `prepare_scripts_for_serverless_exec()` and trusts the result.

---

## 2. Script Generation: prepare_scripts_for_serverless_exec()

**File:** `/home/ubuntu/pash/compiler/serverless/ir_helper.py`

**Location:** Lines 307-370

### Function Signature

```python
def prepare_scripts_for_serverless_exec(
    ir: IR,                          # Optimized IR graph
    shell_vars: dict,                # Shell variables
    args: argparse.Namespace,        # PaSh arguments
    declared_functions_filename,     # Path to declared functions file
    recover: bool = False            # Recovery mode flag
) -> Tuple[str, str, Dict[str, str], Set[str]]:
    """
    Returns:
        main_graph_script_id: Script ID for main shell (EC2 coordinator)
        main_subgraph_script_id: Script ID for first Lambda
        script_id_to_script: Dict mapping script IDs to script content
        ec2_set: Set of script IDs that should run on EC2
    """
```

### What It Does

This is the **core function** that:
1. Splits the IR into subgraphs
2. Augments subgraphs with remote communication nodes
3. Converts each subgraph to a bash script
4. Classifies scripts as EC2 or Lambda

### Code Flow

```python
def prepare_scripts_for_serverless_exec(...):
    # STEP 1: Split IR into subgraphs
    # ================================
    subgraphs, mapping = split_ir(ir)

    # split_ir() returns:
    # - subgraphs: List of IR objects, each representing a computation stage
    # - mapping: Dict mapping edge IDs to subgraphs (for connecting stages)

    # STEP 2: Add remote communication nodes to subgraphs
    # ====================================================
    main_graph_script_id, subgraph_script_id_pairs, main_subgraph_script_id, fifo_to_be_replaced = \
        add_nodes_to_subgraphs(
            subgraphs,              # List of IR subgraphs
            ir.get_file_id_gen(),   # File ID generator (for creating new edges)
            mapping,                # Edge-to-subgraph mapping
            args,                   # PaSh arguments
            recover=recover         # Recovery mode flag
        )

    # add_nodes_to_subgraphs() returns:
    # - main_graph_script_id: UUID for EC2 coordinator script
    # - subgraph_script_id_pairs: Dict mapping each subgraph to a UUID
    # - main_subgraph_script_id: UUID for the first Lambda script
    # - fifo_to_be_replaced: FIFO renaming for recovery mode

    # STEP 3: Load declared functions
    # ================================
    declared_functions = ""
    with open(declared_functions_filename, "r") as f:
        declared_functions = f.read()

    # STEP 4: Convert subgraphs to bash scripts
    # ==========================================
    script_id_to_script = {}
    ec2_set = set()

    for subgraph, id_ in subgraph_script_id_pairs.items():
        # 4a. Create temp directory setup commands
        dir_set = set()
        for edge in subgraph.all_fids():
            if edge.is_ephemeral():
                dir_set.add(os.path.join(config.PASH_TMP_PREFIX, edge.prefix))

        mk_dirs = "mkdir -p "+config.PASH_TMP_PREFIX+" \n"
        for dir in dir_set:
            mk_dirs += "mkdir -p "+dir+" \n"

        # 4b. Add environment setup
        export_path = "export PATH=$PATH:runtime\n"
        export_rust_trace = "export RUST_BACKTRACE=1\n"
        export_lib_path = "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:runtime/lib\n"
        add_version = "version=$2\n"

        # 4c. Convert IR to shell script using to_shell()
        script = export_path + export_lib_path + export_rust_trace + \
                 add_version + mk_dirs + f"{declared_functions}\n" + \
                 to_shell(subgraph, args)

        # 4d. Handle recovery mode FIFO replacements
        if recover and subgraph in fifo_to_be_replaced:
            for new_fifo, recover_fifo in fifo_to_be_replaced[subgraph]:
                script = script.replace(new_fifo, recover_fifo)

        # 4e. Save script
        script_name = os.path.join(config.PASH_TMP_PREFIX, str(id_))
        script_id_to_script[str(id_)] = script
        with open(script_name, "w") as f:
            f.write(script)

        # 4f. Classify as EC2 or Lambda
        if ("split" in script) or ("s3-put" in script) or \
           ("sort -m" in script) or ("merge" in script):
            ec2_set.add(str(id_))

    return str(main_graph_script_id), str(main_subgraph_script_id), \
           script_id_to_script, ec2_set
```

### Key Insight

This function orchestrates three major transformations:
1. **Split IR** → Multiple subgraphs (one per computation stage)
2. **Augment subgraphs** → Add remote communication nodes
3. **Generate scripts** → Convert IR to bash

---

## 3. Graph Splitting: split_ir()

**File:** `/home/ubuntu/pash/compiler/dspash/ir_helper.py`

**Location:** Lines 77-169

### What It Does

Splits an optimized IR into subgraphs, where each subgraph represents a continuous computation segment between splitter and merger nodes.

### Algorithm

```python
def split_ir(graph: IR) -> Tuple[List[IR], Dict[int, IR]]:
    """
    Example: Given this IR:

                  - tr -- grep -
                /                \
    cat - uniq - split            merge - wc
                \                /
                  - tr -- grep -

    Returns 4 subgraphs:
    1. [cat - uniq - split]
    2. [tr - grep]
    3. [tr - grep]
    4. [merge - wc]
    """

    # 1. Initialize BFS queue with source nodes
    source_node_ids = graph.source_nodes()
    subgraphs = []
    queue = deque([(source, IR({}, {})) for source in source_node_ids])

    # 2. Track visited edges and nodes (DAG traversal)
    visited_edges = set(graph.all_input_fids())
    visited_nodes = set()
    input_fifo_map = defaultdict(list)  # edge_id -> [subgraphs using this edge]

    # 3. BFS traversal
    while queue:
        old_node_id, subgraph = queue.popleft()
        input_fids = graph.get_node_input_fids(old_node_id)
        output_fids = graph.get_node_output_fids(old_node_id)

        # 3a. Check if all inputs are ready
        if any(fid not in visited_edges for fid in input_fids):
            if subgraph.source_nodes():
                subgraphs.append(subgraph)
            continue

        # 3b. SPLIT POINT: If this is a merger (multiple inputs), start new subgraph
        if len(input_fids) > 1 and subgraph.source_nodes():
            if subgraph not in subgraphs:
                subgraphs.append(subgraph)
            subgraph = IR({}, {})  # Create new subgraph

        # 3c. Skip if already visited
        if old_node_id in visited_nodes:
            continue
        visited_nodes.add(old_node_id)

        # 3d. Copy node to current subgraph
        node = graph.get_node(old_node_id).copy()
        node_id = node.get_id()

        # 3e. Add input edges to subgraph
        for input_fid in input_fids:
            if input_fid.get_ident() not in subgraph.edges:
                subgraph.add_to_edge(input_fid, node_id)
            else:
                subgraph.set_edge_to(input_fid.get_ident(), node_id)
            input_fifo_map[input_fid.get_ident()].append(subgraph)

        # 3f. Add output edges to subgraph
        for output_fid in output_fids:
            subgraph.add_from_edge(node_id, output_fid)
            visited_edges.add(output_fid)

        # 3g. Add the node itself
        subgraph.add_node(node)

        # 3h. Continue to next nodes
        next_ids = graph.get_next_nodes(old_node_id)

        if len(next_ids) == 1:
            # Single output: continue in same subgraph
            queue.append((next_ids[0], subgraph))
        else:
            # SPLIT POINT: Multiple outputs (splitter node)
            subgraphs.append(subgraph)
            # Create new subgraph for each output
            for next_id in next_ids:
                queue.append((next_id, IR({}, {})))

    return subgraphs, input_fifo_map
```

### Example: sort with width=2

**Input IR:**
```
s3-get → r_split → [sort, sort] → sort -m → s3-put
```

**Output subgraphs:**
1. `[s3-get → r_split]` - Splitter subgraph
2. `[sort]` - Worker 1
3. `[sort]` - Worker 2
4. `[sort -m → s3-put]` - Merger subgraph

---

## 4. Communication Node Injection: add_nodes_to_subgraphs()

**File:** `/home/ubuntu/pash/compiler/serverless/ir_helper.py`

**Location:** Lines 28-304

### What It Does

This is the **most complex function**. It:
1. Connects subgraphs using remote communication nodes (pashlib)
2. Generates unique communication keys (UUIDs) for each connection
3. Adds sender and receiver nodes to subgraphs
4. Tracks which subgraphs need which pashlib arguments
5. Handles recovery mode with ingates/outgates

### High-Level Algorithm

```python
def add_nodes_to_subgraphs(
    subgraphs: List[IR],
    file_id_gen: FileIdGen,
    input_fifo_map: Dict[int, IR],
    args: argparse.Namespace,
    recover: bool = False
):
    # Data structures
    main_graph = IR({}, {})                    # EC2 coordinator graph
    subgraph_script_id_pairs = {}              # subgraph → UUID mapping
    main_subgraph_script_id = None             # First Lambda UUID
    subgraph_stun_lib_args = {}                # subgraph → [pashlib args]
    key_to_sender_receiver = {}                # UUID → (sender_subgraph, receiver_subgraph)
    key_to_data_type = {}                      # UUID → "batch" or "line"

    # PHASE 1: Replace output edges with remote send
    # ===============================================
    for subgraph in subgraphs:
        sink_nodes = subgraph.sink_nodes()
        out_edges = subgraph.get_node_output_fids(sink_nodes[0])

        for out_edge in out_edges:
            # 1a. Create communication key (UUID)
            communication_key = uuid4()

            # 1b. Replace edge with ephemeral edge
            out_edge_id = out_edge.get_ident()
            ephemeral_edge = file_id_gen.next_ephemeral_file_id()
            subgraph.replace_edge(out_edge_id, ephemeral_edge)

            # 1c. Check if this is the last subgraph (no downstream)
            last_subgraph = (out_edge_id not in input_fifo_map)

            if last_subgraph:
                # Send to S3/stdout
                communication_key = "stdout" or filename
                remote_write = make_serverless_remote_pipe(
                    local_fifo_id=ephemeral_edge.get_ident(),
                    is_remote_read=False,
                    remote_key=communication_key,
                    output_edge=stdout,
                    is_tcp=False  # Use S3
                )
                subgraph.add_node(remote_write)
            else:
                # Send to another subgraph via pashlib
                # Format: send*UUID*sender*receiver*fifo_path
                arg = f"send*{communication_key}*0*1*{config.PASH_TMP_PREFIX}{ephemeral_edge}"

                if subgraph not in subgraph_stun_lib_args:
                    subgraph_stun_lib_args[subgraph] = []
                subgraph_stun_lib_args[subgraph].append(arg)

                # Track sender/receiver
                key_to_sender_receiver[str(communication_key)] = [subgraph, None]

                # Determine data type (batch or line)
                out_node = subgraph.get_node(sink_nodes[0])
                if isinstance(out_node, RWrap) or isinstance(out_node, RSplit):
                    key_to_data_type[str(communication_key)] = "batch"
                else:
                    key_to_data_type[str(communication_key)] = "line"

            # 1d. Find matching downstream subgraph
            if out_edge_id in input_fifo_map and out_edge.is_ephemeral():
                matching_subgraph = input_fifo_map[out_edge_id][0]

                # Assign UUID to downstream subgraph if not already assigned
                if matching_subgraph not in subgraph_script_id_pairs:
                    script_identifier = uuid4()
                    subgraph_script_id_pairs[matching_subgraph] = script_identifier

                # Replace edge in downstream subgraph
                new_edge = file_id_gen.next_file_id()
                new_edge.set_resource(out_edge.get_resource())
                matching_subgraph.replace_edge(out_edge.get_ident(), new_edge)

                # Add receiver argument
                arg = f"recv*{communication_key}*1*0*{config.PASH_TMP_PREFIX}{new_edge}"
                if matching_subgraph not in subgraph_stun_lib_args:
                    subgraph_stun_lib_args[matching_subgraph] = []
                subgraph_stun_lib_args[matching_subgraph].append(arg)

                # Update sender/receiver tracking
                key_to_sender_receiver[str(communication_key)][1] = matching_subgraph

    # PHASE 2: Replace non-ephemeral input edges (files, stdin)
    # ==========================================================
    for subgraph in subgraphs:
        if subgraph not in subgraph_script_id_pairs:
            main_subgraph_script_id = uuid4()
            subgraph_script_id_pairs[subgraph] = main_subgraph_script_id

        source_nodes = subgraph.source_nodes()
        for source in source_nodes:
            for in_edge in subgraph.get_node_input_fids(source):
                if in_edge.has_file_resource():
                    filename = in_edge.get_resource().uri

                    # Replace with ephemeral edge + remote read
                    ephemeral_edge = file_id_gen.next_ephemeral_file_id()
                    subgraph.replace_edge(in_edge.get_ident(), ephemeral_edge)

                    remote_read = make_serverless_remote_pipe(
                        local_fifo_id=ephemeral_edge.get_ident(),
                        is_remote_read=True,
                        remote_key=filename,
                        output_edge=None,
                        is_tcp=False  # Read from S3
                    )
                    subgraph.add_node(remote_read)

    # PHASE 3: Recovery mode - Add ingates/outgates
    # ==============================================
    if recover:
        # Classify subgraphs (lambda, splitter, merger, etc.)
        subgraph_types = classify_subgraphs(subgraph_stun_lib_args)

        # Add versioning to keys and FIFOs
        for subgraph, stun_lib_args in subgraph_stun_lib_args.items():
            if subgraph_types[subgraph] == "lambda":
                # Lambda: key → key_v0
                for i in range(len(stun_lib_args)):
                    key = stun_lib_args[i].split("*")[1]
                    stun_lib_args[i] = stun_lib_args[i].replace(key, key+"_v0")
            else:
                # EC2: Add ingates/outgates for recovery
                for i in range(len(stun_lib_args)):
                    key = stun_lib_args[i].split("*")[1]
                    stun_lib_args[i] = stun_lib_args[i].replace(key, key+"_v0")
                    fifo = stun_lib_args[i].split("*")[4]
                    stun_lib_args[i] = stun_lib_args[i].replace(fifo, fifo+"_v0")

                    if "send" in stun_lib_args[i]:
                        # Add outgate
                        receiver = key_to_sender_receiver[key][1]
                        receiver_script_id = subgraph_script_id_pairs[receiver]
                        outgate_args = [fifo, fifo, key, "900", str(receiver_script_id)]
                        outgate = make_serverless_outgate(key_to_data_type[key], outgate_args)
                        subgraph.add_node(outgate)
                    else:
                        # Add ingate
                        sender = key_to_sender_receiver[key][0]
                        sender_script_id = subgraph_script_id_pairs[sender]
                        ingate_args = [fifo, key, fifo, str(sender_script_id)]
                        ingate = make_serverless_ingate(key_to_data_type[key], ingate_args)
                        subgraph.add_node(ingate)

    # PHASE 4: Create pashlib nodes
    # ==============================
    for subgraph, stun_lib_args in subgraph_stun_lib_args.items():
        # Group arguments by peer (avoid conflicts)
        args_lists = [[]]
        peer_lists = [set()]

        for arg in stun_lib_args:
            key = arg.split("*")[1].split("_v")[0] if recover else arg.split("*")[1]
            sender, receiver = key_to_sender_receiver[key]
            peer = receiver if "send" in arg else sender

            # Find a group where this peer isn't already present
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

        # Create pashlib node for each group
        for args_list in args_lists:
            stun_lib = make_serverless_remote_pipe_one_proc(args_list)
            subgraph.add_node(stun_lib)

    # Create main graph ID
    main_graph_script_id = uuid4()
    subgraph_script_id_pairs[main_graph] = main_graph_script_id

    return main_graph_script_id, subgraph_script_id_pairs, \
           main_subgraph_script_id, fifo_to_be_renamed
```

### Example: Communication Setup for width=2 sort

**After split_ir:**
```
Subgraph 1: [s3-get → r_split]
Subgraph 2: [sort]
Subgraph 3: [sort]
Subgraph 4: [sort -m → s3-put]
```

**After add_nodes_to_subgraphs:**
```
Subgraph 1: [s3-get → r_split → pashlib send*uuid1* send*uuid2*]
Subgraph 2: [pashlib recv*uuid1* → sort → pashlib send*uuid3*]
Subgraph 3: [pashlib recv*uuid2* → sort → pashlib send*uuid4*]
Subgraph 4: [pashlib recv*uuid3* recv*uuid4* → sort -m → s3-put]
```

**Communication keys generated:**
- `uuid1`: Subgraph 1 → Subgraph 2
- `uuid2`: Subgraph 1 → Subgraph 3
- `uuid3`: Subgraph 2 → Subgraph 4
- `uuid4`: Subgraph 3 → Subgraph 4

---

## 5. Script Generation: to_shell()

**File:** `/home/ubuntu/pash/compiler/ir_to_ast.py`

**Function:** `to_shell(subgraph: IR, args: argparse.Namespace) -> str`

### What It Does

Converts an IR graph to bash script. This function:
1. Creates FIFO management functions
2. Generates pipeline commands for each node
3. Handles process synchronization
4. Adds cleanup code

### Generated Script Structure

```bash
#!/bin/bash
export PATH=$PATH:runtime
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:runtime/lib
export RUST_BACKTRACE=1

# Helper functions
pash_communicate_daemon() { ... }

# Create temp directories
mkdir -p /tmp/pash_xxx/hash/

# FIFO management
rm_pash_fifos() {
  rm -f /tmp/pash_xxx/hash/#fifo1
  rm -f /tmp/pash_xxx/hash/#fifo2
}

mkfifo_pash_fifos() {
  mkfifo /tmp/pash_xxx/hash/#fifo1
  mkfifo /tmp/pash_xxx/hash/#fifo2
}

rm_pash_fifos
mkfifo_pash_fifos

pids_to_kill=""

# Pipeline commands (background processes)
{ /opt/pashlib send*uuid1*0*1*/tmp/fifo1 recv*uuid2*1*0*/tmp/fifo2 & }
pids_to_kill="${!} ${pids_to_kill}"

{ python3.9 aws/s3-get-object.py "inputs/1M.txt" /tmp/fifo1 & }
pids_to_kill="${!} ${pids_to_kill}"

{ r_split -r /tmp/fifo2 1000000 /tmp/fifo3 /tmp/fifo4 & }
pids_to_kill="${!} ${pids_to_kill}"

# Wait for completion
source runtime/wait_for_output_and_sigpipe_rest.sh ${pids_to_kill}

rm_pash_fifos
( exit "${internal_exec_status}" )
```

---

## 6. Putting It All Together: Complete Function Call Chain

```
ServerlessManager.run()
    ↓ (receives request via socket)
ServerlessManager.handler(request, conn)
    ↓
init(ir_filename)
    ├─→ read_graph(filename)
    │       └─→ pickle.load(ir_file)  # Returns (ir, shell_vars, args)
    └─→ config.set_config_globals_from_pash_args(args)
    ↓
prepare_scripts_for_serverless_exec(ir, shell_vars, args, declared_functions_filename, recover)
    ├─→ split_ir(ir)  [from dspash/ir_helper.py]
    │       └─→ BFS traversal, split at splitter/merger nodes
    │       └─→ Returns: (subgraphs, input_fifo_map)
    │
    ├─→ add_nodes_to_subgraphs(subgraphs, file_id_gen, input_fifo_map, args, recover)
    │       ├─→ Phase 1: Replace output edges with remote send
    │       │       ├─→ uuid4()  # Generate communication key
    │       │       ├─→ make_serverless_remote_pipe()  # For S3 reads/writes
    │       │       └─→ Build pashlib args: "send*UUID*port*count*fifo"
    │       │
    │       ├─→ Phase 2: Replace input file edges with remote read
    │       │       └─→ make_serverless_remote_pipe()  # For S3 reads
    │       │
    │       ├─→ Phase 3: Recovery mode (if enabled)
    │       │       ├─→ classify_subgraphs()
    │       │       ├─→ make_serverless_ingate()
    │       │       └─→ make_serverless_outgate()
    │       │
    │       └─→ Phase 4: Create pashlib nodes
    │               └─→ make_serverless_remote_pipe_one_proc(args_list)
    │       │
    │       └─→ Returns: (main_graph_script_id, subgraph_script_id_pairs,
    │                     main_subgraph_script_id, fifo_to_be_replaced)
    │
    └─→ For each subgraph:
            ├─→ to_shell(subgraph, args)  [from ir_to_ast.py]
            │       ├─→ Generate FIFO management
            │       ├─→ Convert nodes to commands
            │       └─→ Add synchronization
            │
            ├─→ Save script to /tmp/{uuid}
            └─→ Classify as EC2 or Lambda (based on keywords)
    ↓
    └─→ Returns: (main_graph_script_id, main_subgraph_script_id,
                  script_id_to_script, ec2_set)
    ↓
ServerlessManager.handler (continued)
    ├─→ s3.put_object() for each script  [multi-threaded]
    │       └─→ Upload to s3://bucket/sls-scripts/{timestamp}/{uuid}.sh
    │
    └─→ For each script:
            ├─→ If script_id in ec2_set:
            │       └─→ invoke_ec2() or run_local()
            │               └─→ subprocess.run(["/bin/bash", script_path])
            │
            └─→ Else:
                    └─→ invoke_lambda()
                            └─→ lambda_client.invoke(
                                    FunctionName="lambda",
                                    Payload={"ids": [id], "folder_ids": [folder_id]}
                                )
```

---

## 7. Key Data Structures

### IR (Intermediate Representation)

```python
class IR:
    nodes: Dict[int, Node]  # node_id → Node object
    edges: Dict[int, FileId]  # edge_id → FileId object

    def source_nodes(self) -> List[int]
    def sink_nodes(self) -> List[int]
    def get_node_input_fids(self, node_id: int) -> List[FileId]
    def get_node_output_fids(self, node_id: int) -> List[FileId]
```

### FileId (Edge in dataflow graph)

```python
class FileId:
    ident: int
    resource: Resource  # FileResource, FDResource, EphemeralResource

    def is_ephemeral(self) -> bool
    def has_file_resource(self) -> bool
    def get_resource(self) -> Resource
```

### Node (Computation unit)

```python
class Node:
    id: int
    type: str  # "r_split", "sort", "grep", etc.

    # Examples:
    RSplit: Round-robin splitter
    RWrap: Batch wrapper
    ServerlessRemotePipe: S3 read/write or pashlib TCP
```

---

## 8. Example: Complete Flow for sort width=2

### Initial IR (after optimization)

```
Node 1: s3-get-object (input: file "1M.txt", output: ephemeral_1)
Node 2: r_split (input: ephemeral_1, output: ephemeral_2, ephemeral_3)
Node 3: sort (input: ephemeral_2, output: ephemeral_4)
Node 4: sort (input: ephemeral_3, output: ephemeral_5)
Node 5: sort -m (input: ephemeral_4, ephemeral_5, output: ephemeral_6)
Node 6: s3-put-object (input: ephemeral_6, output: file "result.txt")
```

### After split_ir()

**Subgraph 1:**
```
Node 1: s3-get-object
Node 2: r_split
```

**Subgraph 2:**
```
Node 3: sort
```

**Subgraph 3:**
```
Node 4: sort
```

**Subgraph 4:**
```
Node 5: sort -m
Node 6: s3-put-object
```

### After add_nodes_to_subgraphs()

**Subgraph 1 (UUID: aa11):**
```
Node 1: s3-get-object (output: fifo1)
Node 2: r_split (input: fifo1, output: fifo2, fifo3)
Node 3: pashlib send*bb22*0*1*fifo2 send*cc33*0*1*fifo3
```

**Subgraph 2 (UUID: bb22):**
```
Node 1: pashlib recv*bb22*1*0*fifo4
Node 2: sort (input: fifo4, output: fifo5)
Node 3: pashlib send*dd44*0*1*fifo5
```

**Subgraph 3 (UUID: cc33):**
```
Node 1: pashlib recv*cc33*1*0*fifo6
Node 2: sort (input: fifo6, output: fifo7)
Node 3: pashlib send*ee55*0*1*fifo7
```

**Subgraph 4 (UUID: ff66):**
```
Node 1: pashlib recv*dd44*1*0*fifo8 recv*ee55*1*0*fifo9
Node 2: sort -m (input: fifo8, fifo9, output: fifo10)
Node 3: s3-put-object (input: fifo10)
```

### After to_shell()

**Script aa11 (Splitter - EC2):**
```bash
#!/bin/bash
export PATH=$PATH:runtime
mkfifo /tmp/fifo1 /tmp/fifo2 /tmp/fifo3
{ /opt/pashlib send*bb22*0*1*/tmp/fifo2 send*cc33*0*1*/tmp/fifo3 & }
{ python3.9 aws/s3-get-object.py "inputs/1M.txt" /tmp/fifo1 & }
{ r_split -r /tmp/fifo1 1000000 /tmp/fifo2 /tmp/fifo3 & }
wait
```

**Script bb22 (Worker 1 - Lambda):**
```bash
#!/bin/bash
mkfifo /tmp/fifo4 /tmp/fifo5
{ /opt/pashlib recv*bb22*1*0*/tmp/fifo4 send*dd44*0*1*/tmp/fifo5 & }
{ sort /tmp/fifo4 > /tmp/fifo5 & }
wait
```

**Script cc33 (Worker 2 - Lambda):**
```bash
#!/bin/bash
mkfifo /tmp/fifo6 /tmp/fifo7
{ /opt/pashlib recv*cc33*1*0*/tmp/fifo6 send*ee55*0*1*/tmp/fifo7 & }
{ sort /tmp/fifo6 > /tmp/fifo7 & }
wait
```

**Script ff66 (Merger - EC2):**
```bash
#!/bin/bash
mkfifo /tmp/fifo8 /tmp/fifo9 /tmp/fifo10
{ /opt/pashlib recv*dd44*1*0*/tmp/fifo8 recv*ee55*1*0*/tmp/fifo9 & }
{ sort -m /tmp/fifo8 /tmp/fifo9 > /tmp/fifo10 & }
{ python3.9 aws/s3-put-object.py "outputs/result.txt" /tmp/fifo10 & }
wait
```

### Execution

1. **ServerlessManager uploads scripts to S3:**
   - `s3://bucket/sls-scripts/1732578912/aa11.sh`
   - `s3://bucket/sls-scripts/1732578912/bb22.sh`
   - `s3://bucket/sls-scripts/1732578912/cc33.sh`
   - `s3://bucket/sls-scripts/1732578912/ff66.sh`

2. **ServerlessManager invokes workers:**
   - `run_local(aa11)` - EC2 splitter
   - `invoke_lambda([1732578912], [bb22])` - Lambda worker 1
   - `invoke_lambda([1732578912], [cc33])` - Lambda worker 2
   - `run_local(ff66)` - EC2 merger

3. **Workers execute:**
   - EC2 runs `aa11.sh` locally
   - Lambda 1 downloads `bb22.sh` from S3, executes
   - Lambda 2 downloads `cc33.sh` from S3, executes
   - EC2 runs `ff66.sh` locally

4. **Communication:**
   - pashlib establishes TCP connections using UUIDs as rendezvous keys
   - Data streams: EC2 → Lambda 1, EC2 → Lambda 2, Lambda 1 → EC2, Lambda 2 → EC2

---

## 9. Summary: Key Functions and Their Roles

| Function | File | Purpose |
|----------|------|---------|
| `ServerlessManager.handler()` | serverless_executor.py | Orchestrates entire serverless execution |
| `prepare_scripts_for_serverless_exec()` | serverless/ir_helper.py | Converts IR to bash scripts |
| `split_ir()` | dspash/ir_helper.py | Splits IR into subgraphs at split/merge points |
| `add_nodes_to_subgraphs()` | serverless/ir_helper.py | Adds remote communication nodes to subgraphs |
| `to_shell()` | ir_to_ast.py | Converts IR graph to bash script |
| `invoke_lambda()` | serverless_executor.py | Invokes AWS Lambda with script ID |
| `run_local()` | serverless_executor.py | Executes script locally on EC2 |

---

## 10. Critical Insights

### Why split_ir() splits at mergers

Mergers have multiple inputs, meaning multiple upstream computations must complete before the merger can start. This is a natural synchronization point for distributed execution.

### Why add_nodes_to_subgraphs() is complex

It must:
1. Connect subgraphs that were split apart
2. Handle both TCP (pashlib) and S3 communication
3. Generate unique UUIDs for peer-to-peer connections
4. Avoid port conflicts (group by peer)
5. Support recovery mode with versioned FIFOs

### Why EC2 runs split/merge operations

- **Split:** Needs to read entire input from S3, better on EC2 with persistent connection
- **Merge:** Receives streams from multiple Lambdas, easier to coordinate on EC2
- **Lambda:** Best for embarrassingly parallel operations (sort, grep individual chunks)

---

**For new sessions:** Read this document to understand the internal workings of PaSh serverless execution, from IR to deployed Lambda/EC2 scripts.
