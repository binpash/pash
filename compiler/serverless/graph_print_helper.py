
from typing import Dict, List, Tuple
from ir import *
import config

def _get_node_label(node, max_width=60) -> str:
    """
    Extract a readable label for a node including command name and key details.

    Args:
        node: DFGNode instance
        max_width: Maximum width for the label (default 60 chars)

    Returns:
        Human-readable string representation of the node
    """
    cmd_inv = node.cmd_invocation_with_io_vars
    cmd_name = str(cmd_inv.cmd_name)

    # Get basename for cleaner display
    cmd_basename = os.path.basename(cmd_name)

    # Special handling for s3-shared-read - show full arguments without truncation
    if 's3-shard-read' in cmd_basename or 's3-shar' in cmd_basename:
        import re
        parts = [cmd_basename]

        # Extract filename and byte range from operands
        if hasattr(cmd_inv, 'operand_list'):
            for op in cmd_inv.operand_list:
                op_str = str(op)
                # Look for filename (files typically have extensions or paths)
                if any(ext in op_str for ext in ['.csv', '.txt', '.json', 'inputs/', 'outputs/']):
                    filename = op_str.strip('"').strip("'")
                    if filename and not filename.startswith('bytes='):
                        parts.append(f'"{filename}"')
                # Look for byte range
                elif 'bytes=' in op_str:
                    match = re.search(r'bytes=(\d+-\d+)', op_str)
                    if match:
                        parts.append(match.group(0))
                # Look for shard index
                elif 'shard=' in op_str:
                    match = re.search(r'shard=(\d+)', op_str)
                    if match:
                        parts.append(match.group(0))
                # Look for num_shards
                elif 'num_shards=' in op_str:
                    match = re.search(r'num_shards=(\d+)', op_str)
                    if match:
                        parts.append(match.group(0))

        return ' '.join(parts)  # Return without truncation

    # Special handling for r_wrap - extract inner command
    if cmd_basename == "r_wrap" or "wrap" in cmd_basename.lower():
        # Try to extract the inner command from operands
        if hasattr(cmd_inv, 'operand_list') and cmd_inv.operand_list:
            for op in cmd_inv.operand_list:
                op_str = str(op)
                # Look for bash -c followed by command
                if "bash" in op_str and "-c" in op_str:
                    # Skip this one, look at next
                    continue
                # Check if this looks like a command (has spaces or quotes)
                if not op_str.isdigit() and ('"' in op_str or "'" in op_str or " " in op_str):
                    # Extract command, remove quotes
                    inner_cmd = op_str.strip().strip('"').strip("'")
                    # Get first command name
                    first_word = inner_cmd.split()[0] if inner_cmd.split() else inner_cmd
                    # Truncate if needed
                    if len(inner_cmd) > 40:
                        inner_cmd = inner_cmd[:37] + "..."
                    return f"wrap: {first_word}"
        return "wrap"

    # Add flags if present
    flags = []
    if hasattr(cmd_inv, 'flag_option_list') and cmd_inv.flag_option_list:
        flags = [str(flag.get_name()) if hasattr(flag, 'get_name') else str(flag)
                 for flag in cmd_inv.flag_option_list[:3]]  # Limit to first 3 flags

    # Add key operands (excluding edge IDs which are usually integers)
    operands = []
    if hasattr(cmd_inv, 'operand_list') and cmd_inv.operand_list:
        for op in cmd_inv.operand_list[:2]:  # Limit to first 2 operands
            op_str = str(op)
            # Skip if it looks like an edge ID (pure integer)
            if not op_str.isdigit():
                # Detect class repr and simplify
                if "<class '" in op_str or "pash_annotations" in op_str:
                    # Extract just the class name
                    if "'" in op_str:
                        parts = op_str.split("'")
                        if len(parts) > 1:
                            class_path = parts[1]
                            op_str = class_path.split(".")[-1]  # Get last part
                operands.append(op_str if len(op_str) <= 30 else op_str[:27] + "...")

    # Build label
    label = cmd_basename
    if flags:
        label += " " + " ".join(flags)
    if operands:
        label += " " + " ".join(operands)

    # Apply width cap with smart truncation
    if len(label) > max_width:
        # Keep start and end, truncate middle
        keep_chars = (max_width - 5) // 2  # -5 for " ... "
        label = label[:keep_chars] + " ... " + label[-keep_chars:]

    return label


def _get_edge_info(fid, edge_id, subgraph, is_input):
    """
    Extract comprehensive information about an edge for display.

    Args:
        fid: FileId object for this edge
        edge_id: The edge identifier
        subgraph: IR subgraph containing this edge
        is_input: True if this is an input edge, False if output

    Returns:
        Formatted string like:
        - "edge_16 [FILE] \"covid-mts/inputs/in_tiny.csv\""
        - "edge_17 [FIFO] '/tmp/pash_fifo_123'"
        - "edge_18 [stdin]"
    """
    info_parts = [f"edge_{edge_id}"]

    # Edge type indicator
    if fid.is_ephemeral():
        info_parts.append("[FIFO]")
    elif fid.has_file_descriptor_resource():
        resource = fid.get_resource()
        if resource.is_stdin():
            info_parts.append("[stdin]")
        elif resource.is_stdout():
            info_parts.append("[stdout]")
        else:
            info_parts.append("[fd]")
    elif fid.has_file_resource():
        info_parts.append("[FILE]")

    # Resource details
    if fid.has_file_resource():
        filename = str(fid.get_resource().uri).strip('"')
        info_parts.append(f'"{filename}"')

        # Check if this edge connects to an S3 node with byte range info
        # Look for ServerlessRemotePipe nodes in this subgraph
        for node_id in subgraph.nodes.keys():
            node = subgraph.get_node(node_id)
            # Check if this node reads from this edge
            if edge_id in node.get_input_list():
                # Check if it's an S3 lambda node
                if hasattr(node, 'cmd_invocation_with_io_vars'):
                    cmd = node.cmd_invocation_with_io_vars
                    # Look for byte range in operands
                    if hasattr(cmd, 'operand_list'):
                        for op in cmd.operand_list:
                            op_str = str(op)
                            if 'bytes=' in op_str:
                                # Extract byte range
                                import re
                                match = re.search(r'bytes=(\d+-\d+)', op_str)
                                if match:
                                    info_parts.append(match.group(0))
                                    break  # Found byte range, no need to continue
    elif fid.is_ephemeral():
        # Show FIFO name if needed
        fifo_suffix = fid.get_fifo_suffix()
        if fifo_suffix:
            info_parts.append(f"'{config.PASH_TMP_PREFIX}{fifo_suffix}'")

    return " ".join(info_parts)


def _find_edge_connections(subgraphs: List[IR]) -> Dict[int, Tuple[int, int]]:
    """
    Build a map of cross-subgraph edge connections.

    Args:
        subgraphs: List of IR subgraphs

    Returns:
        Dictionary mapping edge_id -> (source_subgraph_idx, dest_subgraph_idx)
        where -1 indicates external input/output
    """
    edge_to_sg = {}  # edge_id -> subgraph_idx
    edge_connections = {}  # edge_id -> (from_sg_idx, to_sg_idx)

    # First pass: map each edge to its subgraph
    for sg_idx, subgraph in enumerate(subgraphs):
        for edge_id in subgraph.edges.keys():
            if edge_id not in edge_to_sg:
                edge_to_sg[edge_id] = []
            edge_to_sg[edge_id].append(sg_idx)

    # Second pass: identify cross-subgraph edges
    for edge_id, sg_indices in edge_to_sg.items():
        if len(sg_indices) > 1:
            # Edge connects multiple subgraphs
            # Find which subgraph outputs to this edge and which inputs from it
            for sg_idx, subgraph in enumerate(subgraphs):
                if edge_id in subgraph.edges:
                    _, from_node, to_node = subgraph.edges[edge_id]

                    # If this subgraph has a node outputting to this edge
                    if from_node is not None and to_node is None:
                        # This is a source subgraph for this edge
                        for other_idx in sg_indices:
                            if other_idx != sg_idx:
                                edge_connections[edge_id] = (sg_idx, other_idx)

    return edge_connections


def _draw_subgraph(sg_idx: int, subgraph: IR, edge_connections: Dict[int, Tuple[int, int]]) -> List[str]:
    """
    Generate ASCII art for a single subgraph.

    Args:
        sg_idx: Index of this subgraph
        subgraph: IR object representing the subgraph
        edge_connections: Map of edge_id -> (from_sg, to_sg)

    Returns:
        List of strings (lines) representing the ASCII visualization
    """
    lines = []

    # Header
    node_count = len(subgraph.nodes)
    lines.append(f"Subgraph {sg_idx} ({node_count} node{'s' if node_count != 1 else ''}):")

    if node_count == 0:
        lines.append("  (empty)")
        return lines

    # Get source nodes (nodes with no incoming edges from other nodes in subgraph)
    source_nodes = subgraph.source_nodes()

    # Build a simple topological view
    visited = set()

    def draw_node_chain(node_id, indent=2):
        """Recursively draw nodes in execution order"""
        if node_id in visited:
            return
        visited.add(node_id)

        node = subgraph.get_node(node_id)
        label = _get_node_label(node)

        # Get input edges
        input_ids = node.get_input_list()
        output_ids = node.get_output_list()

        # Check if this node has multiple external inputs (merge pattern)
        external_inputs = []
        for in_id in input_ids:
            if in_id in subgraph.edges:
                fid, from_node, _ = subgraph.edges[in_id]
                if from_node is None:  # External input
                    external_inputs.append((in_id, fid))

        # Show converging inputs for merge nodes
        if len(external_inputs) > 1:
            # Draw merge pattern: inputs converge into the node
            box_width = max(len(label) + 4, 12)
            for i, (in_id, fid) in enumerate(external_inputs):
                source_info = ""
                if in_id in edge_connections:
                    from_sg, _ = edge_connections[in_id]
                    source_info = f" (from Subgraph {from_sg})"
                elif fid.has_file_resource():
                    source_info = f" (file: {fid.get_resource().uri})"

                if i == 0:
                    # First input
                    lines.append(" " * indent + f"[edge_{in_id}{source_info}] ─┐")
                elif i == len(external_inputs) - 1:
                    # Last input
                    lines.append(" " * (indent + box_width + 3) + "├─> +" + "-" * (box_width - 2) + "+")
                    lines.append(" " * indent + f"[edge_{in_id}{source_info}] ─┘   | " + label.ljust(box_width - 4) + " |")
                    lines.append(" " * (indent + box_width + 7) + "+" + "-" * (box_width - 2) + "+")
                else:
                    # Middle inputs
                    lines.append(" " * (indent + box_width + 3) + "│")
                    lines.append(" " * indent + f"[edge_{in_id}{source_info}] ─┤")
        else:
            # Single or no external input - show normally
            for in_id in input_ids:
                if in_id in subgraph.edges:
                    fid, from_node, _ = subgraph.edges[in_id]
                    if from_node is None:  # External input
                        source_info = ""
                        if in_id in edge_connections:
                            from_sg, _ = edge_connections[in_id]
                            source_info = f" (from Subgraph {from_sg})"
                        elif fid.has_file_resource():
                            source_info = f" (file: {fid.get_resource().uri})"
                        lines.append(" " * indent + f"[Input: edge_{in_id}{source_info}]")
                        lines.append(" " * indent + "    |")
                        lines.append(" " * indent + "    v")

            # Draw the node box
            box_width = max(len(label) + 4, 12)
            lines.append(" " * indent + "+" + "-" * (box_width - 2) + "+")
            lines.append(" " * indent + "| " + label.ljust(box_width - 4) + " |")
            lines.append(" " * indent + "+" + "-" * (box_width - 2) + "+")

        # Get next nodes
        next_nodes = subgraph.get_next_nodes(node_id)

        # Check if this node has multiple outputs (split pattern)
        has_multiple_outputs = len(output_ids) > 1

        if has_multiple_outputs:
            # Draw split pattern with branches (always show as branches if multiple outputs)
            lines.append(" " * indent + "    |")
            for i, out_id in enumerate(output_ids):
                dest_info = ""
                if out_id in edge_connections:
                    _, to_sg = edge_connections[out_id]
                    dest_info = f" (to Subgraph {to_sg})"
                elif out_id in subgraph.edges:
                    fid, _, to_node = subgraph.edges[out_id]
                    if to_node is None and fid.has_file_resource():
                        dest_info = f" (file: {fid.get_resource().uri})"

                branch_marker = "├──>" if i < len(output_ids) - 1 else "└──>"
                lines.append(" " * indent + f"    {branch_marker} [edge_{out_id}{dest_info}]")
        elif len(next_nodes) == 0:
            # Sink node with single output
            for out_id in output_ids:
                lines.append(" " * indent + "    |")
                lines.append(" " * indent + "    v")
                dest_info = ""
                if out_id in edge_connections:
                    _, to_sg = edge_connections[out_id]
                    dest_info = f" (to Subgraph {to_sg})"
                elif out_id in subgraph.edges:
                    fid, _, to_node = subgraph.edges[out_id]
                    if to_node is None and fid.has_file_resource():
                        dest_info = f" (file: {fid.get_resource().uri})"
                lines.append(" " * indent + f"[Output: edge_{out_id}{dest_info}]")
        elif len(next_nodes) == 1:
            # Linear chain continues
            lines.append(" " * indent + "    |")
            lines.append(" " * indent + "    v")
            draw_node_chain(next_nodes[0], indent)
        else:
            # Multiple next nodes (shouldn't happen if has_multiple_outputs handled above)
            lines.append(" " * indent + "    |")
            for i, next_id in enumerate(next_nodes):
                edge_label = ""
                for out_id in output_ids:
                    if out_id in subgraph.edges:
                        _, _, to_node = subgraph.edges[out_id]
                        if to_node == next_id:
                            edge_label = f"edge_{out_id}"
                            if out_id in edge_connections:
                                _, to_sg = edge_connections[out_id]
                                edge_label += f" (to Subgraph {to_sg})"
                            break

                branch_marker = "├──>" if i < len(next_nodes) - 1 else "└──>"
                lines.append(" " * indent + f"    {branch_marker} [{edge_label}]")

    # Start from source nodes
    for source_id in source_nodes:
        draw_node_chain(source_id)
        lines.append("")  # Blank line between chains

    return lines


def _create_subgraph_summary(sg_idx: int, subgraph: IR) -> str:
    """
    Create a one-line summary of a subgraph for compact display.

    Args:
        sg_idx: Index of the subgraph
        subgraph: IR object

    Returns:
        String summary like "cat → r_split → [2 outputs]"
    """
    if len(subgraph.nodes) == 0:
        return "(empty)"

    # NEW APPROACH: Show all nodes starting with source nodes
    # (Can't follow edges for s3 nodes since they have output_edge=None)

    # Start with source nodes (nodes that start the pipeline)
    source_nodes = subgraph.source_nodes()
    node_ids_in_order = list(source_nodes)  # Source nodes first

    # Add remaining nodes in sorted order
    for node_id in sorted(subgraph.nodes.keys()):
        if node_id not in node_ids_in_order:
            node_ids_in_order.append(node_id)

    # Create labels
    node_labels = []
    for node_id in node_ids_in_order:
        node = subgraph.get_node(node_id)
        label = _get_node_label(node, max_width=100)
        node_labels.append(label)
        if len(node_labels) >= 5:  # Limit to 5 nodes
            break

    if len(node_labels) == 0:
        return "(no nodes)"

    # Check if there are multiple outputs at the end
    sink_nodes = subgraph.sink_nodes()
    if len(sink_nodes) > 0:
        total_outputs = sum(len(subgraph.get_node_output_fids(s)) for s in sink_nodes)
        if total_outputs > 1:
            node_labels.append(f"[{total_outputs} outputs]")

    return " → ".join(node_labels)


def _draw_unified_graph(subgraphs: List[IR], edge_connections: Dict[int, Tuple[int, int]]) -> List[str]:
    """
    Draw all subgraphs in a unified view with bridge connections.

    Args:
        subgraphs: List of IR subgraphs
        edge_connections: Map of edge_id -> (from_sg, to_sg)

    Returns:
        List of lines representing the unified visualization
    """
    lines = []

    # Build a mapping of which subgraphs connect to which
    sg_connections = {}  # sg_idx -> {to_sg_idx: [edge_ids]}
    for edge_id, (from_sg, to_sg) in edge_connections.items():
        if from_sg not in sg_connections:
            sg_connections[from_sg] = {}
        if to_sg not in sg_connections[from_sg]:
            sg_connections[from_sg][to_sg] = []
        sg_connections[from_sg][to_sg].append(edge_id)

    # Draw each subgraph with connections
    for sg_idx, subgraph in enumerate(subgraphs):
        node_count = len(subgraph.nodes)
        summary = _create_subgraph_summary(sg_idx, subgraph)

        # NEW: Find and display external INPUT edges
        source_nodes = subgraph.source_nodes()
        for source_id in source_nodes:
            input_fids = subgraph.get_node_input_fids(source_id)
            for fid in input_fids:
                edge_id = fid.get_ident()
                if edge_id in subgraph.edges:
                    _, from_node, _ = subgraph.edges[edge_id]
                    if from_node is None:  # External input
                        edge_info = _get_edge_info(fid, edge_id, subgraph, is_input=True)
                        lines.append(f"[Input: {edge_info}]")
                        lines.append("    ↓")

        # Create box around subgraph
        title = f" Subgraph {sg_idx} ({node_count} node{'s' if node_count != 1 else ''}) "
        box_width = max(len(title) + 4, len(summary) + 6, 50)

        # Top border
        lines.append("┌" + "─" * (len(title)) + "─" + "─" * (box_width - len(title) - 2) + "┐")
        lines.append("│" + title + " " * (box_width - len(title) - 2) + "│")
        lines.append("│" + " " * (box_width - 2) + "│")

        # Content
        lines.append("│  " + summary + " " * (box_width - len(summary) - 4) + "│")
        lines.append("│" + " " * (box_width - 2) + "│")

        # Bottom border
        lines.append("└" + "─" * (box_width - 2) + "┘")

        # NEW: Find and display external OUTPUT edges
        sink_nodes = subgraph.sink_nodes()
        has_external_output = False
        for sink_id in sink_nodes:
            output_fids = subgraph.get_node_output_fids(sink_id)
            for fid in output_fids:
                edge_id = fid.get_ident()
                if edge_id in subgraph.edges:
                    _, _, to_node = subgraph.edges[edge_id]
                    if to_node is None and edge_id not in edge_connections:
                        # External output (not to another subgraph)
                        if not has_external_output:
                            lines.append("    ↓")
                            has_external_output = True
                        edge_info = _get_edge_info(fid, edge_id, subgraph, is_input=False)
                        lines.append(f"[Output: {edge_info}]")

        # Draw bridge connections to next subgraphs
        if sg_idx in sg_connections:
            for to_sg_idx, edge_ids in sg_connections[sg_idx].items():
                # Draw dashed line(s) to next subgraph
                for edge_id in edge_ids:
                    connector = f"    ╎ edge_{edge_id} ─ ─ ─ ─> (to Subgraph {to_sg_idx})"
                    lines.append(connector)
            lines.append("")  # Blank line after connections
        else:
            lines.append("")  # Blank line between subgraphs

    return lines



def pretty_print_subgraphs(subgraphs: List[IR], show_connections: bool = True, unified_view: bool = False) -> None:
    """
    Print a human-readable ASCII graph representation of subgraphs.

    This function visualizes the structure of IR subgraphs with ASCII art,
    showing nodes as boxes, edges as arrows, and cross-subgraph connections.

    Args:
        subgraphs: List of IR objects (subgraphs from split_ir)
        show_connections: Whether to show edge connections summary at end (separate view only)
        unified_view: If True, show all subgraphs in a connected layout with bridge edges

    Example:
        >>> subgraphs, mapping = split_ir(optimized_ir)
        >>> pretty_print_subgraphs(subgraphs)  # Detailed separate view
        >>> pretty_print_subgraphs(subgraphs, unified_view=True)  # Compact unified view
    """
    # Find all cross-subgraph connections
    edge_connections = _find_edge_connections(subgraphs)

    if unified_view:
        # Unified view: all subgraphs in connected layout
        print("=" * 80)
        print(" " * 20 + "UNIFIED GRAPH VISUALIZATION")
        print("=" * 80)
        print()

        lines = _draw_unified_graph(subgraphs, edge_connections)
        for line in lines:
            print(line)
    else:
        # Separate view: each subgraph shown individually
        print("=" * 80)
        print(" " * 25 + "SUBGRAPH VISUALIZATION")
        print("=" * 80)
        print()

        # Draw each subgraph
        for sg_idx, subgraph in enumerate(subgraphs):
            lines = _draw_subgraph(sg_idx, subgraph, edge_connections)
            for line in lines:
                print(line)
            print()

        # Show cross-subgraph connection summary
        if show_connections and edge_connections:
            print("=" * 80)
            print(" " * 20 + "CROSS-SUBGRAPH CONNECTIONS")
            print("=" * 80)
            for edge_id, (from_sg, to_sg) in sorted(edge_connections.items()):
                print(f"  edge_{edge_id}: Subgraph {from_sg} → Subgraph {to_sg}")
            print()
