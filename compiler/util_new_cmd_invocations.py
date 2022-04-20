# imports from annotation framework
import sys
sys.path.insert(1, "/home/felix/git-repos/MIT/annotations")
# for typing
from datatypes_new.CommandInvocationPrefix import CommandInvocationPrefix
# for use
from definitions.ir.dfg_node import DFGNode

def get_command_invocation_prefix_from_dfg_node(dfg_node):
    return CommandInvocationPrefix(cmd_name = dfg_node.com_name,
                                   flag_option_list = dfg_node.flag_option_list,
                                   positional_config_list = dfg_node.positional_config_list)

