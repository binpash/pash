# splitter, merger: change key to key_v0, fifo in the pashlib to fifo_v0, for each pashlib, know the downstream script id, add command
# lambda: change key to key_v${version}

import argparse
from copy import deepcopy
from typing import Dict, List, Tuple
import sys
import os
from uuid import uuid4
sys.path.append(os.path.join(os.getenv("PASH_TOP"), "compiler"))
import definitions.ir.nodes.serverless_remote_pipe as serverless_remote_pipe
import definitions.ir.nodes.serverless_lambda_invoke as serverless_lambda_invoke
from definitions.ir.nodes.r_wrap import RWrap
from definitions.ir.nodes.r_split import RSplit
from dspash.ir_helper import split_ir
from ir_to_ast import to_shell
from ir import *
import config
import pash_compiler


# Debugpy setup - connect before main logic runs
debug = False
import debugpy
if debug:
    debugpy.listen(5683)
    print("🐛 Debugger listening on port 5683. Attach VSCode now!")
    debugpy.wait_for_client()
    print("✅ Debugger attached! Continuing execution...")

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
    # Convert to set as we only need to do this once per unique edge/lambda
    rsplit_output_id_set = set(rsplit_output_ids)

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
                    # YES - Replace with S3 direct read
                    subgraph.replace_edge(in_edge_id, ec2_file_in_edge)

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
    

# cat -> split
def optimize_s3_lambda_direct_streaming(subgraphs:List[IR]):
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
    assert len(source_nodes) == 1
    #TODO. encode assertions to run on narrow execution path 
    for source in source_nodes:
        in_edges = first_subgraph.get_node_input_fids(source)
        assert len(in_edges) == 1
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
                        return subgraphs
                    next_node_id = next_nodes_ids[0]
                    next_node = first_subgraph.get_node(next_node_id)

                    if isinstance(next_node, RSplit):
                        # Extract rsplit output edge IDs
                        #can do in a more straightforward way with a helper method in ir
                        rsplit_output_fids = first_subgraph.get_node_output_fids(next_node_id)
                        rsplit_output_ids = [fid.get_ident() for fid in rsplit_output_fids]

                        # Call adjust function with extracted info
                        # TODO. dont need first subgraph anymore
                        adjust_lambda_incoming_edges(
                            first_subgraph,
                            subgraphs[1:],
                            file_id_gen,
                            s3_file_path,       # S3 file path
                            rsplit_output_ids,  # rsplit output edge IDs,
                            in_edge
                        )

                        return subgraphs[1:]

    return subgraphs

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
    key_to_data_type = {} # batch or line
    eager_edges = []

    #modify the first subgraph if it is a s3 -> cat -> split
    # then the loop will not even have it anymore
    # helper fn
    
    #subgraphs = optimize_s3_lambda_direct_streaming(subgraphs)


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
                    subgraph.add_node(remote_read) # TODO. investigate if this makes a node with a file or s3 input?

                    pash_compiler.add_eager(ephemeral_edge.get_ident(), subgraph, file_id_gen)
                
                else:
                    # sometimes a command can have both a file resource and an ephemeral resources (example: spell oneliner)
                    continue
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
        add_version = "version=$2\n"
        script = export_path+export_lib_path+export_rust_trace+add_version+mk_dirs+f"{declared_functions}\n"+to_shell(subgraph, args)
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
