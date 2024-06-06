import argparse
from typing import Dict, List, Tuple
import sys
import os
from uuid import uuid4
sys.path.append(os.path.join(os.getenv("PASH_TOP"), "compiler"))
import definitions.ir.nodes.serverless_remote_pipe as serverless_remote_pipe
import definitions.ir.nodes.serverless_lambda_invoke as serverless_lambda_invoke
from definitions.ir.nodes.r_merge import RMerge
from dspash.ir_helper import split_ir
from ir_to_ast import to_shell
from ir import *
import config


def add_stdout_fid(graph : IR, file_id_gen: FileIdGen) -> FileId:
    stdout = file_id_gen.next_file_id()
    stdout.set_resource(FileDescriptorResource(('fd', 1)))
    graph.add_edge(stdout)
    return stdout

def add_nodes_to_subgraphs(subgraphs:List[IR], file_id_gen: FileIdGen, input_fifo_map:Dict[int, IR], args: argparse.Namespace):
    """ Takes a list of subgraphs and augments subgraphs with the necessary remote
        read/write nodes for data movement and lambda invocation nodes to trigger
        downstream processing. This function also produces graph that should run in
        the original shell in which pash was executed. This graph contains
        remote read/write nodes for stdin/stdout, named pipes, and files.

    Args:
        subgraphs: list of sub sections of an optimized IR (returned from split_ir)
        file_id_gen: file id generator of the original ir
        input_fifo_map: mapping from input idge id to subgraph (returned from split_ir)
    Returns:
        main_graph_script_id: the script id to execute on main shell
        subgraph_script_id_pairs: mapping from subgraph to unique script id
        main_subgraph_script_id: the script id to execute in the first lambda
    """
    # The graph to execute in the main pash_compiler
    main_graph = IR({}, {})
    subgraph_script_id_pairs = {}
    main_subgraph_script_id = None

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
                    subgraph.add_node(serverless_lambda_invoke.make_serverless_lambda_invoke(script_identifier))
            else:
                # Add edge to main graph
                matching_subgraph = main_graph
                matching_subgraph.add_edge(new_edge)
            remote_read = serverless_remote_pipe.make_serverless_remote_pipe(local_fifo_id=new_edge.get_ident(),
                                                                              is_remote_read=True,
                                                                              remote_key=communication_key,
                                                                              output_edge=new_edge,
                                                                              is_tcp=(not matching_subgraph is main_graph))
            matching_subgraph.add_node(remote_read)

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
                    remote_read = serverless_remote_pipe.make_serverless_remote_pipe(local_fifo_id=ephemeral_edge.get_ident(),
                                                                              is_remote_read=True,
                                                                              remote_key=filename,
                                                                              output_edge=None,
                                                                              is_tcp=False)
                    subgraph.add_node(remote_read)
                else:
                    # sometimes a command can have both a file resource and an ephemeral resources (example: spell oneliner)
                    continue
    main_graph_script_id = uuid4()
    subgraph_script_id_pairs[main_graph] = main_graph_script_id
    return main_graph_script_id, subgraph_script_id_pairs, main_subgraph_script_id


def prepare_scripts_for_serverless_exec(ir: IR, shell_vars: dict, args: argparse.Namespace) -> Tuple[str, str, Dict[str, str]]:
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
    main_graph_script_id, subgraph_script_id_pairs, main_subgraph_script_id = add_nodes_to_subgraphs(subgraphs, ir.get_file_id_gen(), mapping, args)

    # save the output scripts
    script_id_to_script = {}
    for subgraph, id_ in subgraph_script_id_pairs.items():
        # making necessary temp directories
        dir_set = set()
        for edge in subgraph.all_fids():
            if edge.is_ephemeral():
                dir_set.add(os.path.join(config.PASH_TMP_PREFIX, edge.prefix))
        mk_dirs = "mkdir -p "+config.PASH_TMP_PREFIX+" \n"
        for dir in dir_set:
            mk_dirs += "mkdir -p "+dir+" \n"
        script = mk_dirs+to_shell(subgraph, args)
        # generate scripts
        script_name = os.path.join(config.PASH_TMP_PREFIX, str(id_))
        script_id_to_script[str(id_)] = script
        with open (script_name, "w") as f:
            f.write(script)
        if id_ == main_graph_script_id:
            log("Script for main shell saved in:"+script_name)
        elif id_ == main_subgraph_script_id:
            log("Script for first lambda saved in:"+script_name)
        else:
            log("Script for other lambda saved in:"+script_name)

    return str(main_graph_script_id), str(main_subgraph_script_id), script_id_to_script
