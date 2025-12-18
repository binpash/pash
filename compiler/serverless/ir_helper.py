# splitter, merger: change key to key_v0, fifo in the pashlib to fifo_v0, for each pashlib, know the downstream script id, add command
# lambda: change key to key_v${version}

import argparse
from copy import deepcopy
from typing import Dict, List, Tuple
import sys
import os
from uuid import uuid4

import boto3
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

# cat -> split
def optimize_s3_lambda_direct_streaming(subgraphs:List[IR]):
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

                        return subgraphs[1:], in_edge, total_lambdas

    return subgraphs, None, 0

def change_sorting(subgraphs:List[IR]) -> List[IR]:
    for subgraph in subgraphs:
        source_nodes = subgraph.source_nodes()    # list of ints
        for source in source_nodes:
            source_node = subgraph.get_node(source)
            if source_node.cmd_name == 'sort': #todo also do this with sort -m
                source_node.cmd_name = 'LC_ALL=C sort'


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
    Returns:
        main_graph_script_id: the script id to execute on main shell
        subgraph_script_id_pairs: mapping from subgraph to unique script id
        main_subgraph_script_id: the script id to execute in the first lambda
    """
    # Print subgraphs for debugging/visualization
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
    subgraphs, ec2_in_edge, total_lambdas = optimize_s3_lambda_direct_streaming(subgraphs)
    print("AFTER OPTIMIZATION, TOTAL LAMBDAS:", total_lambdas)


    print("\n" + "="*80)
    print("DEBUG: Visualizing subgraphs after S3 optimization")
    print("="*80)

    pretty_print_subgraphs(subgraphs, unified_view=True)

    print("="*80)
    exit()

    
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

    lambda_counter = 0
    job_uid = uuid4()
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



                        # # Get basename and strip extension and quotes
                        # basename = str(filename).split("/")[-1].replace('.txt"', "").replace('"', '')

                        # # Parse size with suffix (e.g., "1M" -> 1048576)
                        # filesize = None
                        # for suffix, multiplier in size_multipliers.items():
                        #     if basename.endswith(suffix):
                        #         filesize = int(basename.replace(suffix, "")) * multiplier
                        #         break

                        # # Fallback: try to parse as plain integer
                        # if filesize is None:
                        #     try:
                        #         filesize = int(basename)
                        #     except ValueError:
                        #         # Default to 1MB if parsing fails
                        #         print(f"Warning: Could not parse file size from '{filename}', defaulting to 1MB")
                        #         assert False, "File size parsing failed" # important to know if this is possibly failing, also byte range will give incorrect results here

                        print(f"File: {filename} -> Size: {filesize} bytes")

                        # Calculate byte range for this lambda
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
                                                                                job_uid=job_uid)
                        
                        # s3getobj byterange=start-end -> sort   : lambda1
                        lambda_counter += 1
                    else:
                        remote_read = serverless_remote_pipe.make_serverless_remote_pipe(local_fifo_id=ephemeral_edge.get_ident(),
                                                                                is_remote_read=True,
                                                                                remote_key=filename,
                                                                                output_edge=None,
                                                                                is_tcp=False)
                    subgraph.add_node(remote_read) # This makes a node with an s3 input (converts the input filename to an s3 get obj cmd)
                
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

    # Get nodes in execution order
    source_nodes = subgraph.source_nodes()
    if not source_nodes:
        return "(no source)"

    # Build a simple chain of command names
    visited = set()
    chain = []

    def build_chain(node_id):
        if node_id in visited or len(chain) >= 5:  # Limit to 5 nodes
            return
        visited.add(node_id)

        node = subgraph.get_node(node_id)
        label = _get_node_label(node, max_width=20)
        chain.append(label)

        next_nodes = subgraph.get_next_nodes(node_id)
        if len(next_nodes) == 1:
            build_chain(next_nodes[0])
        elif len(next_nodes) > 1:
            chain.append(f"[{len(next_nodes)} outputs]")

    build_chain(source_nodes[0])

    if len(chain) == 0:
        return "(no nodes)"

    return " → ".join(chain)


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
