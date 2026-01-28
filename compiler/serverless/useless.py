def add_tee_nodes_after_rsplit(
    first_subgraph: IR,
    file_id_gen: FileIdGen,
    s3_output_prefix: str = "outputs/"
) -> bool:
    """
    Adds tee nodes after rsplit in the first subgraph to save copies to S3.

    This function:
    1. Finds the rsplit node in the first subgraph
    2. For each rsplit output:
       - Creates intermediate temp edge
       - Inserts tee node to duplicate the stream
       - Creates s3-put-object node to save copy to S3
    3. Rewires the graph to maintain original dataflow

    Args:
        first_subgraph: The first subgraph (typically contains rsplit)
        file_id_gen: FileIdGen to create new edge IDs
        s3_output_prefix: S3 prefix for output files (default "outputs/")

    Returns:
        True if tee nodes were added, False if rsplit not found or other error
    """
    print("\n" + "="*80)
    print("[DEBUG add_tee_nodes_after_rsplit] STARTING")
    print("="*80)

    # Find rsplit node in the first subgraph
    rsplit_node = None
    rsplit_node_id = None

    print(f"[DEBUG] Searching for RSplit node in first_subgraph")
    print(f"[DEBUG] first_subgraph has {len(first_subgraph.nodes)} nodes")

    for node_id in first_subgraph.nodes.keys():
        node = first_subgraph.get_node(node_id)
        print(f"[DEBUG] Checking node {node_id}: {type(node).__name__}")
        if isinstance(node, RSplit):
            rsplit_node = node
            rsplit_node_id = node_id
            print(f"[DEBUG] Found RSplit node with id: {rsplit_node_id}")
            break

    if rsplit_node is None:
        print("[DEBUG] No rsplit found, returning False")
        return False

    # Get rsplit output edges
    rsplit_output_fids = first_subgraph.get_node_output_fids(rsplit_node_id)
    print(f"[DEBUG] RSplit has {len(rsplit_output_fids)} output edges")

    if len(rsplit_output_fids) == 0:
        print("[DEBUG] No outputs, returning False")
        return False

    # For each rsplit output edge, insert tee + s3-put
    for idx, rsplit_output_fid in enumerate(rsplit_output_fids):
        print(f"\n[DEBUG] Processing rsplit output {idx+1}/{len(rsplit_output_fids)}")
        original_edge_id = rsplit_output_fid.get_ident()
        print(f"[DEBUG]   original_edge_id: {original_edge_id}")

        # Get the current consumer of the original edge
        _fid, from_node, to_node = first_subgraph.edges[original_edge_id]
        print(f"[DEBUG]   Edge info: from_node={from_node}, to_node={to_node}")

        # 1. Create temp edge (between rsplit and tee)
        print(f"[DEBUG] Step 1: Creating temp edge")
        temp_edge = file_id_gen.next_ephemeral_file_id()
        print(f"[DEBUG]   temp_edge id: {temp_edge.get_ident()}")
        first_subgraph.add_edge(temp_edge)

        # 2. Create s3 upload edge (between tee and s3-put)
        print(f"[DEBUG] Step 2: Creating s3 edge")
        s3_edge = file_id_gen.next_ephemeral_file_id()
        print(f"[DEBUG]   s3_edge id: {s3_edge.get_ident()}")
        first_subgraph.add_edge(s3_edge)

        # 3. Rewire rsplit → temp_edge (instead of rsplit → original_edge)
        print(f"[DEBUG] Step 3: Rewiring rsplit → temp_edge")
        rsplit_node.replace_edge(original_edge_id, temp_edge.get_ident())
        first_subgraph.set_edge_from(rsplit_node_id, temp_edge)

        # 4. Create and add tee node
        print(f"[DEBUG] Step 4: Creating tee node")
        tee_node = make_tee_node(
            temp_edge.get_ident(),      # Input: from rsplit
            original_edge_id,            # Output 1: to downstream (original path)
            s3_edge.get_ident()         # Output 2: to s3-put
        )
        print(f"[DEBUG]   Adding tee node to subgraph")
        first_subgraph.add_node(tee_node)

        # 5. Connect edges to/from tee
        print(f"[DEBUG] Step 5: Connecting edges to/from tee")
        first_subgraph.set_edge_to(temp_edge.get_ident(), tee_node.get_id())       # temp → tee (input)

        # Get the actual FileId from the graph's edges for the original edge
        original_fid = first_subgraph.edges[original_edge_id][0]
        print(f"[DEBUG]   Retrieved original_fid from graph: {original_fid}")

        first_subgraph.set_edge_from(original_fid, tee_node.get_id())               # tee → original (stdout)
        print("HEREHERE")
        first_subgraph.set_edge_from(s3_edge, tee_node.get_id())                    # tee → s3 (file arg)
        print(f"[DEBUG]   Edges connected")

        # 6. Original edge now flows from tee (not rsplit)
        # Keep the downstream connection intact: original_edge → to_node
        print(f"[DEBUG] Step 6: Setting original edge downstream connection")
        first_subgraph.set_edge_to(original_edge_id, to_node)

        # 7. Create and add s3-put node
        print(f"[DEBUG] Step 7: Creating s3-put node")
        fifo_name = f"fifo{original_edge_id}"
        s3_key = f"{s3_output_prefix}{fifo_name}.out"
        print(f"[DEBUG]   s3_key: {s3_key}")

        s3_put_node = make_s3_put_node(
            s3_edge.get_ident(),
            s3_key
        )
        print(f"[DEBUG]   Adding s3-put node to subgraph")
        first_subgraph.add_node(s3_put_node)
        first_subgraph.set_edge_to(s3_edge.get_ident(), s3_put_node.get_id())
        print(f"[DEBUG]   s3-put node connected")

    print("\n" + "="*80)
    print("[DEBUG add_tee_nodes_after_rsplit] COMPLETED SUCCESSFULLY")
    print("="*80 + "\n")
    return True


def make_tee_node(input_id, original_output_id, s3_output_id):
    """
    Creates a tee node that reads from input and writes to both outputs.

    This duplicates the stream:
    - Reads from input_id (rsplit temp output)
    - Writes to original_output_id via stdout (for downstream subgraph)
    - Writes to s3_output_id as file argument (for s3-put-object)

    Args:
        input_id: Input fifo edge ID
        original_output_id: Output edge ID for downstream processing
        s3_output_id: Output edge ID for file/S3 upload

    Returns:
        DFGNode configured as tee command
    """
    print(f"[DEBUG make_tee_node] Creating tee node:")
    print(f"  input_id: {input_id} (type: {type(input_id)})")
    print(f"  original_output_id: {original_output_id} (type: {type(original_output_id)})")
    print(f"  s3_output_id: {s3_output_id} (type: {type(s3_output_id)})")

    access_map = {
        input_id: make_stream_input(),           # Input from rsplit
        original_output_id: make_stream_output(), # Stdout → original fifo
        s3_output_id: make_stream_output()        # File arg → s3 upload fifo
    }
    print(f"[DEBUG make_tee_node] access_map created: {access_map}")

    # tee command: tee <file_operand>
    # Reads stdin, writes to stdout AND to <file_operand>
    operand_list = [s3_output_id]  # File argument for S3 copy
    print(f"[DEBUG make_tee_node] operand_list: {operand_list}")

    cmd_inv = CommandInvocationWithIOVars(
        cmd_name="tee",
        flag_option_list=[],
        operand_list=operand_list,
        implicit_use_of_streaming_input=input_id,
        implicit_use_of_streaming_output=original_output_id,
        access_map=access_map
    )
    print(f"[DEBUG make_tee_node] cmd_inv created: {cmd_inv}")

    node = DFGNode(cmd_inv)
    print(f"[DEBUG make_tee_node] DFGNode created with id: {node.get_id()}")
    return node




def make_s3_put_node(input_id, s3_key, version_arg="$1"):
    """
    Creates an s3-put-object node.

    Command: python3.9 aws/s3-put-object.py <s3_key> <input_fifo> <version>

    Args:
        input_id: Input fifo edge ID
        s3_key: S3 key for the uploaded object
        version_arg: Version argument passed to s3-put-object.py (default "$1")

    Returns:
        DFGNode configured for aws/s3-put-object.py
    """
    print(f"[DEBUG make_s3_put_node] Creating s3-put node:")
    print(f"  input_id: {input_id} (type: {type(input_id)})")
    print(f"  s3_key: {s3_key}")
    print(f"  version_arg: {version_arg}")

    access_map = {input_id: make_stream_input()}
    print(f"[DEBUG make_s3_put_node] access_map: {access_map}")

    operand_list = [
        Operand(Arg.string_to_arg(s3_key)),      # S3 key
        input_id,                                 # Input fifo
        Operand(Arg.string_to_arg(version_arg))  # Version ($1)
    ]
    print(f"[DEBUG make_s3_put_node] operand_list: {operand_list}")

    full_operand_list = [Operand(Arg.string_to_arg("aws/s3-put-object.py"))] + operand_list
    print(f"[DEBUG make_s3_put_node] full_operand_list: {full_operand_list}")

    cmd_inv = CommandInvocationWithIOVars(
        cmd_name="python3.9",
        flag_option_list=[],
        operand_list=full_operand_list,
        implicit_use_of_streaming_input=None,
        implicit_use_of_streaming_output=None,
        access_map=access_map
    )
    print(f"[DEBUG make_s3_put_node] cmd_inv created: {cmd_inv}")

    node = DFGNode(cmd_inv)
    print(f"[DEBUG make_s3_put_node] DFGNode created with id: {node.get_id()}")
    return node



def change_sorting(subgraphs:List[IR]) -> List[IR]:
    for subgraph in subgraphs:
        source_nodes = subgraph.source_nodes()    # list of ints
        for source in source_nodes:
            source_node = subgraph.get_node(source)
            if source_node.cmd_name == 'sort': #todo also do this with sort -m
                source_node.cmd_name = 'LC_ALL=C sort'
