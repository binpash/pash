#!/usr/bin/env python3
"""
Test script for pretty_print_subgraphs function.

This script demonstrates how to use the pretty_print_subgraphs function
to visualize IR subgraphs in a human-readable ASCII format.
"""

import sys
import os

# Add compiler to path
sys.path.append(os.path.join(os.path.dirname(__file__), "compiler"))

from compiler.serverless.ir_helper import pretty_print_subgraphs
from compiler.dspash.ir_helper import split_ir
from compiler.ir import IR
from compiler.definitions.ir.file_id import FileIdGen, FileId
from compiler.definitions.ir.resource import FileResource, EphemeralResource
from compiler.definitions.ir.dfg_node import DFGNode
from pash_annotations.datatypes.CommandInvocationWithIOVars import CommandInvocationWithIOVars

def create_test_ir():
    """
    Create a simple test IR representing: cat file.txt | sed 's/foo/bar/' | sort
    """
    file_id_gen = FileIdGen()

    # Create file IDs
    input_file = file_id_gen.next_file_id()
    input_file.set_resource(FileResource("test_input.txt"))

    edge_1 = file_id_gen.next_ephemeral_file_id()
    edge_2 = file_id_gen.next_ephemeral_file_id()

    output_file = file_id_gen.next_file_id()
    output_file.set_resource(FileResource("test_output.txt"))

    # Create nodes
    # Node 1: cat
    cat_cmd = CommandInvocationWithIOVars(
        cmd_name="cat",
        flag_option_list=[],
        operand_list=[input_file.get_ident()],
        implicit_use_of_streaming_input=None,
        implicit_use_of_streaming_output=edge_1.get_ident()
    )
    cat_node = DFGNode(cat_cmd)

    # Node 2: sed
    sed_cmd = CommandInvocationWithIOVars(
        cmd_name="sed",
        flag_option_list=[],
        operand_list=["s/foo/bar/"],
        implicit_use_of_streaming_input=edge_1.get_ident(),
        implicit_use_of_streaming_output=edge_2.get_ident()
    )
    sed_node = DFGNode(sed_cmd)

    # Node 3: sort
    sort_cmd = CommandInvocationWithIOVars(
        cmd_name="sort",
        flag_option_list=[],
        operand_list=[],
        implicit_use_of_streaming_input=edge_2.get_ident(),
        implicit_use_of_streaming_output=output_file.get_ident()
    )
    sort_node = DFGNode(sort_cmd)

    # Build IR
    nodes = {
        cat_node.get_id(): cat_node,
        sed_node.get_id(): sed_node,
        sort_node.get_id(): sort_node
    }

    edges = {
        input_file.get_ident(): (input_file, None, cat_node.get_id()),
        edge_1.get_ident(): (edge_1, cat_node.get_id(), sed_node.get_id()),
        edge_2.get_ident(): (edge_2, sed_node.get_id(), sort_node.get_id()),
        output_file.get_ident(): (output_file, sort_node.get_id(), None)
    }

    ir = IR(nodes, edges, background=False)
    return ir

def main():
    print("Creating test IR...")
    test_ir = create_test_ir()

    print("\nSplitting IR into subgraphs...")
    subgraphs, mapping = split_ir(test_ir)

    print(f"\nGenerated {len(subgraphs)} subgraph(s)\n")

    print("=" * 80)
    print("Testing pretty_print_subgraphs function:")
    print("=" * 80)
    print()

    # Call the pretty print function
    pretty_print_subgraphs(subgraphs, show_connections=True)

    print("\n✓ Test completed successfully!")

if __name__ == "__main__":
    main()
