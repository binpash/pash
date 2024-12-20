import argparse
from copy import deepcopy
from typing import Dict, List, Tuple
import sys
import os
from uuid import uuid4
sys.path.append(os.path.join(os.getenv("PASH_TOP"), "compiler"))
import definitions.ir.nodes.serverless_remote_pipe as serverless_remote_pipe
import definitions.ir.nodes.serverless_lambda_invoke as serverless_lambda_invoke
from dspash.ir_helper import split_ir
from ir_to_ast import to_shell
from ir import *
import config
import pash_compiler


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

    subgraph_to_downstream_subgraphs = {} # subgraph -> [all downstream subgraphs], used for level order traversal
    subgraph_to_invocation_node = {} # subgraph -> the invocation node invoking this subgraph
    subgraphs_not_having_upstream = set(subgraphs)
    subgraph_stun_lib_args = {}
    key_to_sender_receiver = {}

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
                    
                if str(communication_key) not in key_to_sender_receiver:
                    key_to_sender_receiver[str(communication_key)] = [subgraph, None]
                else:
                    key_to_sender_receiver[str(communication_key)][1] = subgraph

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
                matching_subgraph = main_graph
                matching_subgraph.add_edge(new_edge)
            if (not matching_subgraph is main_graph):
                # arg = recv_rdvkey_1_0_output-fifoname
                arg = "recv*"+str(communication_key)+"*1*0*"+config.PASH_TMP_PREFIX+str(new_edge)
                if matching_subgraph not in subgraph_stun_lib_args:
                    subgraph_stun_lib_args[matching_subgraph] = []
                subgraph_stun_lib_args[matching_subgraph].append(arg)

                if str(communication_key) not in key_to_sender_receiver:
                    key_to_sender_receiver[str(communication_key)] = [None, matching_subgraph]
                else:
                    key_to_sender_receiver[str(communication_key)][1] = matching_subgraph

                if not args.no_eager:
                    pash_compiler.add_eager(new_edge.get_ident(), matching_subgraph, file_id_gen)

            else:
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

    for subgraph, stun_lib_args in subgraph_stun_lib_args.items():
        # stun_lib = serverless_remote_pipe.make_serverless_remote_pipe_one_proc(stun_lib_args)
        # subgraph.add_node(stun_lib)

        args_lists = [[]]
        peer_lists = [set()]
        for arg in stun_lib_args:
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

    if args.sls_instance == "hybrid":
        log("[Serverless Manager] Hybrid instance selected, running non-parallized subgraphs on ec2")
        # level order traversal
        cur_level = list(subgraphs_not_having_upstream)
        next_level = []
        levels = []
        visited = deepcopy(subgraphs_not_having_upstream)
        levels.append([item for item in cur_level])
        while cur_level:
            subgraph = cur_level.pop(0)
            for downstream_subgraph in subgraph_to_downstream_subgraphs.get(subgraph, []):
                if downstream_subgraph not in visited:
                    next_level.append(downstream_subgraph)
                    visited.add(downstream_subgraph)
            if not cur_level:
                cur_level = next_level
                next_level = []
                levels.append([item for item in cur_level])
        for level in levels:
            if len(level) == 1:
                # only one subgraph in this level, need to make it run on ec2
                subgraph = level[0]
                if subgraph in subgraph_to_invocation_node:
                    incovation_node = subgraph_to_invocation_node[subgraph]
                    incovation_node.change_instance("ec2")

    return main_graph_script_id, subgraph_script_id_pairs, main_subgraph_script_id


def prepare_scripts_for_serverless_exec(ir: IR, shell_vars: dict, args: argparse.Namespace, declared_functions_filename) -> Tuple[str, str, Dict[str, str]]:
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
        script = export_path+export_lib_path+export_rust_trace+mk_dirs+f"{declared_functions}\n"+to_shell(subgraph, args)
        # generate scripts
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
