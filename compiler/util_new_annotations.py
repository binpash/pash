# imports from annotation framework
import sys
sys.path.insert(1, "/home/felix/git-repos/MIT/annotations")
# for typing
from datatypes_new.CommandInvocation import CommandInvocation
from annotation_generation_new.datatypes.InputOutputInfo import InputOutputInfo
from annotation_generation_new.datatypes.ParallelizabilityInfo import ParallelizabilityInfo
# for use
from annotation_generation_new.AnnotationGeneration import get_input_output_info_from_cmd_invocation, \
    get_parallelizability_info_from_cmd_invocation

def get_input_output_info_from_cmd_invocation_util(cmd_invocation : CommandInvocation) -> InputOutputInfo:
    return get_input_output_info_from_cmd_invocation(cmd_invocation)

def get_parallelizability_info_from_cmd_invocation_util(cmd_invocation : CommandInvocation) -> ParallelizabilityInfo:
    return get_parallelizability_info_from_cmd_invocation(cmd_invocation)

