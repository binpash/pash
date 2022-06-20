# imports from annotation framework
import sys
sys.path.insert(1, "/home/felix/git-repos/MIT/annotations")
# for typing
from datatypes_new.CommandInvocationInitial import CommandInvocationInitial
from annotation_generation_new.datatypes.InputOutputInfo import InputOutputInfo
from annotation_generation_new.datatypes.ParallelizabilityInfo import ParallelizabilityInfo
from annotation_generation_new.datatypes.CommandProperties import CommandProperties
# for use
from annotation_generation_new.AnnotationGeneration import get_input_output_info_from_cmd_invocation, \
    get_parallelizability_info_from_cmd_invocation

def get_input_output_info_from_cmd_invocation_util(cmd_invocationInitial : CommandInvocationInitial) -> InputOutputInfo:
    return get_input_output_info_from_cmd_invocation(cmd_invocationInitial)

def get_parallelizability_info_from_cmd_invocation_util(cmd_invocationInitial : CommandInvocationInitial) -> ParallelizabilityInfo:
    return get_parallelizability_info_from_cmd_invocation(cmd_invocationInitial)

def construct_property_container_from_list_of_properties(list_properties):
    return CommandProperties(dict(list_properties))

