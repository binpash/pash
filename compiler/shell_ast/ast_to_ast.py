from enum import Enum
import pickle

import config

from shell_ast.ast_util import *
from parse import from_ast_objects_to_shell
from speculative import util_spec

## There are two types of ast_to_ast transformations
class TransformationType(Enum):
    PASH = 'pash'
    SPECULATIVE = 'spec'

## Use this object to pass options inside the preprocessing
## trasnformation.
class TransformationOptions:
    def __init__(self, mode: TransformationType):
        self.mode = mode
            
    def get_mode(self):
        return self.mode


## TODO: Turn it into a Transformation State class, and make a subclass for
##       each of the two transformations. It is important for it to be state, because
##       it will need to be passed around while traversing the tree.
class SpeculativeTransformationState(TransformationOptions):
    def __init__(self, mode: TransformationType, po_file: str):
        super().__init__(mode)
        assert(self.mode is TransformationType.SPECULATIVE)
        self.partial_order_file = po_file

    def get_partial_order_file(self):
        assert(self.mode is TransformationType.SPECULATIVE)
        return self.partial_order_file



##
## Preprocessing
##

## The preprocessing pass replaces all _candidate_ dataflow regions with
## calls to PaSh's runtime to let it establish if they are actually dataflow
## regions. The pass serializes all candidate dataflow regions:
## - A list of ASTs if at the top level or
## - an AST subtree if at a lower level
##
## The PaSh runtime then deserializes the(m, compiles them (if safe) and optimizes them.

preprocess_cases = {
    "Pipe": (lambda trans_options, last_object:
             lambda ast_node: preprocess_node_pipe(ast_node, trans_options, last_object=last_object)),
    "Command": (lambda trans_options, last_object:
                lambda ast_node: preprocess_node_command(ast_node, trans_options, last_object=last_object)),
    "Redir": (lambda trans_options, last_object:
              lambda ast_node: preprocess_node_redir(ast_node, trans_options, last_object=last_object)),
    "Background": (lambda trans_options, last_object:
                   lambda ast_node: preprocess_node_background(ast_node, trans_options, last_object=last_object)),
    "Subshell": (lambda trans_options, last_object:
                   lambda ast_node: preprocess_node_subshell(ast_node, trans_options, last_object=last_object)),
    "For": (lambda trans_options, last_object:
            lambda ast_node: preprocess_node_for(ast_node, trans_options, last_object=last_object)),
    "While": (lambda trans_options, last_object:
              lambda ast_node: preprocess_node_while(ast_node, trans_options, last_object=last_object)),
    "Defun": (lambda trans_options, last_object:
              lambda ast_node: preprocess_node_defun(ast_node, trans_options, last_object=last_object)),
    "Semi": (lambda trans_options, last_object:
             lambda ast_node: preprocess_node_semi(ast_node, trans_options, last_object=last_object)),
    "Or": (lambda trans_options, last_object:
           lambda ast_node: preprocess_node_or(ast_node, trans_options, last_object=last_object)),
    "And": (lambda trans_options, last_object:
            lambda ast_node: preprocess_node_and(ast_node, trans_options, last_object=last_object)),
    "Not": (lambda trans_options, last_object:
            lambda ast_node: preprocess_node_not(ast_node, trans_options, last_object=last_object)),
    "If": (lambda trans_options, last_object:
            lambda ast_node: preprocess_node_if(ast_node, trans_options, last_object=last_object)),
    "Case": (lambda trans_options, last_object:
             lambda ast_node: preprocess_node_case(ast_node, trans_options, last_object=last_object))
}



## Replace candidate dataflow AST regions with calls to PaSh's runtime.
def replace_ast_regions(ast_objects, trans_options):

    preprocessed_asts = []
    candidate_dataflow_region = []
    last_object = False
    for i, ast_object in enumerate(ast_objects):
        # log("Preprocessing AST {}".format(i))
        # log(ast_object)
        ## If we are working on the last object we need to keep that in mind when replacing.
        ##
        ## The last df-region should not be executed in parallel no matter what (to not lose its exit code.)
        if (i == len(ast_objects) - 1):
            # log("Last object")
            last_object = True

        ast, original_text, _linno_before, _linno_after = ast_object
        ## TODO: Turn the untyped ast to an AstNode

        ## Goals: This transformation can approximate in several directions.
        ##        1. Not replacing a candidate dataflow region.
        ##        2. Replacing a too large candidate region
        ##           (making expansion not happen as late as possible)
        ##        3. Not replacing a maximal dataflow region,
        ##           e.g. splitting a big one into two.
        ##        4. Replacing sections that are *certainly* not dataflow regions.
        ##           (This can only lead to performance issues.)
        ##
        ##        Which of the above can we hope to be precise with?
        ##        Can we have proofs indicating that we are not approximating those?

        ## Preprocess ast by replacing subtrees with calls to runtime.
        ## - If the whole AST needs to be replaced (e.g. if it is a pipeline)
        ##   then the second output is true.
        ## - If the next AST needs to be replaced too (e.g. if the current one is a background)
        ##   then the third output is true
        preprocessed_ast_object = preprocess_node(ast, trans_options, last_object=last_object)
        ## If the dataflow region is not maximal then it implies that the whole
        ## AST should be replaced.
        assert(not preprocessed_ast_object.is_non_maximal() 
               or preprocessed_ast_object.should_replace_whole_ast())
        
        ## If the whole AST needs to be replaced then it implies that
        ## something will be replaced
        assert(not preprocessed_ast_object.should_replace_whole_ast() 
               or preprocessed_ast_object.will_anything_be_replaced())

        ## If it isn't maximal then we just add it to the candidate
        if(preprocessed_ast_object.is_non_maximal()):
            candidate_dataflow_region.append((preprocessed_ast_object.ast,
                                              original_text))
        else:
            ## If the current candidate dataflow region is non-empty
            ## it means that the previous AST was in the background so
            ## the current one has to be included in the process no matter what
            if (len(candidate_dataflow_region) > 0):
                candidate_dataflow_region.append((preprocessed_ast_object.ast,
                                                  original_text))
                ## Since the current one is maximal (or not wholy replaced)
                ## we close the candidate.
                dataflow_region_asts, dataflow_region_lines = unzip(candidate_dataflow_region)
                dataflow_region_text = join_original_text_lines(dataflow_region_lines)
                replaced_ast = replace_df_region(dataflow_region_asts, trans_options,
                                                 ast_text=dataflow_region_text, disable_parallel_pipelines=last_object)
                candidate_dataflow_region = []
                preprocessed_asts.append(replaced_ast)
            else:
                if(preprocessed_ast_object.should_replace_whole_ast()):
                    replaced_ast = replace_df_region([preprocessed_ast_object.ast], trans_options,
                                                     ast_text=original_text, disable_parallel_pipelines=last_object)
                    preprocessed_asts.append(replaced_ast)
                else:
                    ## In this case, it is possible that no replacement happened,
                    ## meaning that we can simply return the original parsed text as it was.
                    if(preprocessed_ast_object.will_anything_be_replaced() or original_text is None):
                        preprocessed_asts.append(preprocessed_ast_object.ast)
                    else:
                        preprocessed_asts.append(UnparsedScript(original_text))

    ## Close the final dataflow region
    if(len(candidate_dataflow_region) > 0):
        dataflow_region_asts, dataflow_region_lines = unzip(candidate_dataflow_region)
        dataflow_region_text = join_original_text_lines(dataflow_region_lines)
        replaced_ast = replace_df_region(dataflow_region_asts, trans_options,
                                         ast_text=dataflow_region_text, disable_parallel_pipelines=True)
        candidate_dataflow_region = []
        preprocessed_asts.append(replaced_ast)

    return preprocessed_asts

## This function joins original unparsed shell source in a safe way 
##   so as to deal with the case where some of the text is None (e.g., in case of stdin parsing).
def join_original_text_lines(shell_source_lines_or_none):
    if any([text_or_none is None for text_or_none in shell_source_lines_or_none]):
        return None
    else:
        return "\n".join(shell_source_lines_or_none)

def preprocess_node(ast_object, trans_options, last_object=False):
    global preprocess_cases
    return ast_match_untyped(ast_object, preprocess_cases, trans_options, last_object)

## This preprocesses the AST node and also replaces it if it needs replacement .
## It is called by constructs that cannot be included in a dataflow region.
def preprocess_close_node(ast_object, trans_options, last_object=False):
    preprocessed_ast_object = preprocess_node(ast_object, trans_options, last_object=last_object)
    preprocessed_ast = preprocessed_ast_object.ast
    should_replace_whole_ast = preprocessed_ast_object.should_replace_whole_ast()
    if(should_replace_whole_ast):
        final_ast = replace_df_region([preprocessed_ast], trans_options,
                                      disable_parallel_pipelines=last_object)
        something_replaced = True
    else:
        final_ast = preprocessed_ast
        something_replaced = preprocessed_ast_object.will_anything_be_replaced()
    return final_ast, something_replaced

def preprocess_node_pipe(ast_node, trans_options, last_object=False):
    ## A pipeline is *always* a candidate dataflow region.
    ## Q: Is that true?

    ## TODO: Preprocess the internals of the pipe to allow
    ##       for mutually recursive calls to PaSh.
    ##
    ##       For example, if a command in the pipe has a command substitution
    ##       in one of its arguments then we would like to call our runtime
    ##       there instead of
    preprocessed_ast_object = PreprocessedAST(ast_node,
                                              replace_whole=True,
                                              non_maximal=ast_node.is_background,
                                              last_ast=last_object)
    return preprocessed_ast_object

## TODO: Complete this
def preprocess_node_command(ast_node, trans_options, last_object=False):
    ## TODO: Preprocess the internals of the pipe to allow
    ##       for mutually recursive calls to PaSh.
    ##
    ##       For example, if a command in the pipe has a command substitution
    ##       in one of its arguments then we would like to call our runtime
    ##       there instead of

    ## If there are no arguments, the command is just an
    ## assignment (Q: or just redirections?)
    if(len(ast_node.arguments) == 0):
        preprocessed_ast_object = PreprocessedAST(ast_node,
                                                  replace_whole=False,
                                                  non_maximal=False,
                                                  something_replaced=False,
                                                  last_ast=last_object)
        return preprocessed_ast_object

    ## This means we have a command. Commands are always candidate dataflow
    ## regions.
    preprocessed_ast_object = PreprocessedAST(ast_node,
                                              replace_whole=True,
                                              non_maximal=False,
                                              last_ast=last_object)
    return preprocessed_ast_object

# Background of (linno * t * redirection list) 
## TODO: It might be possible to actually not close the inner node but rather apply the redirections on it
def preprocess_node_redir(ast_node, trans_options, last_object=False):
    preprocessed_node, something_replaced = preprocess_close_node(ast_node.node, trans_options,
                                                                  last_object=last_object)
    ## TODO: Could there be a problem with the in-place update
    ast_node.node = preprocessed_node
    preprocessed_ast_object = PreprocessedAST(ast_node,
                                              replace_whole=False,
                                              non_maximal=False,
                                              something_replaced=something_replaced,
                                              last_ast=last_object)
    return preprocessed_ast_object

## TODO: Is that correct? Also, this should probably affect `semi`, `and`, and `or`
def preprocess_node_background(ast_node, trans_options, last_object=False):
    ## A background node is *always* a candidate dataflow region.
    ## Q: Is that true?

    ## TODO: Preprocess the internals of the background to allow
    ##       for mutually recursive calls to PaSh.
    preprocessed_ast_object = PreprocessedAST(ast_node,
                                              replace_whole=True,
                                              non_maximal=True,
                                              last_ast=last_object)
    return preprocessed_ast_object

## TODO: We can actually preprocess the underlying node and then
##       return its characteristics above. However, we would need
##       to add a field in the IR that a node runs in a subshell
##       (which would have implications on how the backend outputs it).
##
##       e.g. a subshell node should also be output as a subshell in the backend.
## FIXME: This might not just be suboptimal, but also wrong.
def preprocess_node_subshell(ast_node, trans_options, last_object=False):
    preprocessed_body, something_replaced = preprocess_close_node(ast_node.body, trans_options,
                                                                  last_object=last_object)
    ## TODO: Could there be a problem with the in-place update
    ast_node.body = preprocessed_body
    preprocessed_ast_object = PreprocessedAST(ast_node,
                                              replace_whole=False,
                                              non_maximal=False,
                                              something_replaced=something_replaced,
                                              last_ast=last_object)
    return preprocessed_ast_object

## TODO: For all of the constructs below, think whether we are being too conservative

## TODO: This is not efficient at all since it calls the PaSh runtime everytime the loop is entered.
##       We have to find a way to improve that.
def preprocess_node_for(ast_node, trans_options, last_object=False):
    preprocessed_body, something_replaced = preprocess_close_node(ast_node.body, trans_options, last_object=last_object)
    ## TODO: Could there be a problem with the in-place update
    ast_node.body = preprocessed_body
    preprocessed_ast_object = PreprocessedAST(ast_node,
                                              replace_whole=False,
                                              non_maximal=False,
                                              something_replaced=something_replaced,
                                              last_ast=last_object)
    return preprocessed_ast_object

def preprocess_node_while(ast_node, trans_options, last_object=False):
    preprocessed_test, sth_replaced_test = preprocess_close_node(ast_node.test, trans_options, last_object=last_object)
    preprocessed_body, sth_replaced_body = preprocess_close_node(ast_node.body, trans_options, last_object=last_object)
    ## TODO: Could there be a problem with the in-place update
    ast_node.test = preprocessed_test
    ast_node.body = preprocessed_body
    something_replaced = sth_replaced_test or sth_replaced_body
    preprocessed_ast_object = PreprocessedAST(ast_node,
                                              replace_whole=False,
                                              non_maximal=False,
                                              something_replaced=something_replaced,
                                              last_ast=last_object)
    return preprocessed_ast_object

## This is the same as the one for `For`
def preprocess_node_defun(ast_node, trans_options, last_object=False):
    ## TODO: For now we don't want to compile function bodies
    # preprocessed_body = preprocess_close_node(ast_node.body)
    ## TODO: Could there be a problem with the in-place update
    # ast_node.body = preprocessed_body
    preprocessed_ast_object = PreprocessedAST(ast_node,
                                              replace_whole=False,
                                              non_maximal=False,
                                              something_replaced=False,
                                              last_ast=last_object)
    return preprocessed_ast_object

## TODO: If the preprocessed is not maximal we actually need to combine it with the one on the right.
def preprocess_node_semi(ast_node, trans_options, last_object=False):
    # preprocessed_left, should_replace_whole_ast, is_non_maximal = preprocess_node(ast_node.left, irFileGen, config)
    ##
    ## TODO: Is it valid that only the right one is considered the last command?
    preprocessed_left, sth_replaced_left = preprocess_close_node(ast_node.left_operand, trans_options, last_object=False)
    preprocessed_right, sth_replaced_right = preprocess_close_node(ast_node.right_operand, trans_options, last_object=last_object)
    ## TODO: Could there be a problem with the in-place update
    ast_node.left_operand = preprocessed_left
    ast_node.right_operand = preprocessed_right
    sth_replaced = sth_replaced_left or sth_replaced_right
    preprocessed_ast_object = PreprocessedAST(ast_node,
                                              replace_whole=False,
                                              non_maximal=False,
                                              something_replaced=sth_replaced,
                                              last_ast=last_object)
    return preprocessed_ast_object

## TODO: Make sure that what is inside an `&&`, `||`, `!` (and others) does not run in parallel_pipelines 
##       since we need its exit code.
def preprocess_node_and(ast_node, trans_options, last_object=False):
    # preprocessed_left, should_replace_whole_ast, is_non_maximal = preprocess_node(ast_node.left, irFileGen, config)
    preprocessed_left, sth_replaced_left = preprocess_close_node(ast_node.left_operand, trans_options, last_object=last_object)
    preprocessed_right, sth_replaced_right = preprocess_close_node(ast_node.right_operand, trans_options, last_object=last_object)
    ## TODO: Could there be a problem with the in-place update
    ast_node.left_operand = preprocessed_left
    ast_node.right_operand = preprocessed_right
    sth_replaced = sth_replaced_left or sth_replaced_right
    preprocessed_ast_object = PreprocessedAST(ast_node,
                                              replace_whole=False,
                                              non_maximal=False,
                                              something_replaced=sth_replaced,
                                              last_ast=last_object)
    return preprocessed_ast_object

def preprocess_node_or(ast_node, trans_options, last_object=False):
    # preprocessed_left, should_replace_whole_ast, is_non_maximal = preprocess_node(ast_node.left, irFileGen, config)
    preprocessed_left, sth_replaced_left = preprocess_close_node(ast_node.left_operand, trans_options, last_object=last_object)
    preprocessed_right, sth_replaced_right = preprocess_close_node(ast_node.right_operand, trans_options, last_object=last_object)
    ## TODO: Could there be a problem with the in-place update
    ast_node.left_operand = preprocessed_left
    ast_node.right_operand = preprocessed_right
    sth_replaced = sth_replaced_left or sth_replaced_right
    preprocessed_ast_object = PreprocessedAST(ast_node,
                                              replace_whole=False,
                                              non_maximal=False,
                                              something_replaced=sth_replaced,
                                              last_ast=last_object)
    return preprocessed_ast_object

def preprocess_node_not(ast_node, trans_options, last_object=False):
    # preprocessed_left, should_replace_whole_ast, is_non_maximal = preprocess_node(ast_node.left)
    preprocessed_body, sth_replaced = preprocess_close_node(ast_node.body, trans_options, last_object=last_object)
    ## TODO: Could there be a problem with the in-place update
    ast_node.body = preprocessed_body
    preprocessed_ast_object = PreprocessedAST(ast_node,
                                              replace_whole=False,
                                              non_maximal=False,
                                              something_replaced=sth_replaced,
                                              last_ast=last_object)
    return preprocessed_ast_object


def preprocess_node_if(ast_node, trans_options, last_object=False):
    # preprocessed_left, should_replace_whole_ast, is_non_maximal = preprocess_node(ast_node.left, irFileGen, config)
    preprocessed_cond, sth_replaced_cond = preprocess_close_node(ast_node.cond, trans_options, last_object=last_object)
    preprocessed_then, sth_replaced_then = preprocess_close_node(ast_node.then_b, trans_options, last_object=last_object)
    preprocessed_else, sth_replaced_else = preprocess_close_node(ast_node.else_b, trans_options, last_object=last_object)
    ## TODO: Could there be a problem with the in-place update
    ast_node.cond = preprocessed_cond
    ast_node.then_b = preprocessed_then
    ast_node.else_b = preprocessed_else
    sth_replaced = sth_replaced_cond or sth_replaced_then or sth_replaced_else
    preprocessed_ast_object = PreprocessedAST(ast_node,
                                              replace_whole=False,
                                              non_maximal=False,
                                              something_replaced=sth_replaced,
                                              last_ast=last_object)
    return preprocessed_ast_object

def preprocess_case(case, trans_options, last_object=False):
    preprocessed_body, sth_replaced = preprocess_close_node(case["cbody"], trans_options, last_object=last_object)
    case["cbody"] = preprocessed_body
    return case, sth_replaced

def preprocess_node_case(ast_node, trans_options, last_object=False):
    preprocessed_cases_replaced = [preprocess_case(case, trans_options, last_object=last_object) for case in ast_node.cases]
    preprocessed_cases, sth_replaced_cases = list(zip(*preprocessed_cases_replaced))
    ## TODO: Could there be a problem with the in-place update
    ast_node.cases = preprocessed_cases
    preprocessed_ast_object = PreprocessedAST(ast_node,
                                              replace_whole=False,
                                              non_maximal=False,
                                              something_replaced=any(sth_replaced_cases),
                                              last_ast=last_object)
    return preprocessed_ast_object


## TODO: I am a little bit confused about how compilation happens.
##       Does it happen bottom up or top down: i.e. when we first encounter an occurence
##       do we recurse in it and then compile from the leaf, or just compile the surface?



## Replaces IR subtrees with a command that calls them (more
## precisely, a command that calls a python script to call them).
##
## Note: The traversal that replace_irs does, is exactly the same as
## the one that is done by compile_node. Both of these functions
## transform nodes of type t to something else.
##
## TODO: For now this just replaces the IRs starting from the ourside
## one first, but it should start from the bottom up to handle
## recursive IRs.

## This function serializes a candidate df_region in a file, and in its place,
## it adds a command that calls our distribution planner with the name of the
## saved file.
##
## If we are need to disable parallel pipelines, e.g., if we are in the context of an if,
## or if we are in the end of a script, then we set a variable.
def replace_df_region(asts, trans_options, disable_parallel_pipelines=False, ast_text=None):
    transformation_mode = trans_options.get_mode()
    if transformation_mode is TransformationType.PASH:
        ir_filename = ptempfile()

        ## Serialize the node in a file
        with open(ir_filename, "wb") as ir_file:
            pickle.dump(asts, ir_file)

        ## Serialize the candidate df_region asts back to shell
        ## so that the sequential script can be run in parallel to the compilation.
        sequential_script_file_name = ptempfile()
        text_to_output = get_shell_from_ast(asts, ast_text=ast_text)
        ## However, if we have the original ast text, then we can simply output that.
        with open(sequential_script_file_name, "w") as script_file:
            script_file.write(text_to_output)
        replaced_node = make_call_to_pash_runtime(ir_filename, sequential_script_file_name, disable_parallel_pipelines)
    elif transformation_mode is TransformationType.SPECULATIVE:
        ## TODO: This currently writes each command on its own line,
        ##       though it should be improved to better serialize each command in its own file
        ##       and then only saving the ids of each command in the partial order file.
        text_to_output = get_shell_from_ast(asts, ast_text=ast_text)
        ## Generate an ID
        df_region_id = util_spec.get_next_id()
        ## Determine its predecessors
        ## TODO: To make this properly work, we should keep some state
        ##       in the AST traversal to be able to determine predecessors.
        if df_region_id == 0:
            predecessors = []
        else:
            predecessors = [df_region_id - 1]
        ## Write to a file indexed by its ID
        util_spec.save_df_region(text_to_output, trans_options, df_region_id, predecessors)
        ## TODO: Add an entry point to spec through normal PaSh
        replaced_node = make_call_to_spec_runtime(df_region_id)
    else:
        ## Unreachable
        assert(False)

    return replaced_node


def get_shell_from_ast(asts, ast_text=None) -> str:
    ## If we don't have the original ast text, we need to unparse the ast
    if (ast_text is None):
        kv_asts = [ast_node_to_untyped_deep(ast) for ast in asts]
        text_to_output = from_ast_objects_to_shell(kv_asts)
    else:
        text_to_output = ast_text
    return text_to_output


##
## Code that constructs the preprocessed ASTs
##

## TODO: Replace this with the save state
def make_pre_runtime_nodes():
    input_args_command = make_input_args_command()
    save_shell_state_command = make_save_shell_state_command()
    return [save_shell_state_command]

def make_post_runtime_nodes():
    ## TODO: Remove this once pash_runtime is sourced with the same arguments
    set_args_node = restore_arguments_command()
    set_exit_status_node = restore_exit_code_node()
    return [set_exit_status_node]

def make_input_args_command():
    ## Save the input arguments
    ## ```
    ## source $PASH_TOP/runtime/save_args.sh "${@}"
    ## ```
    arguments = [string_to_argument("source"),
                 string_to_argument(config.SAVE_ARGS_EXECUTABLE),
                 [make_quoted_variable("@")]]
    input_args_command = make_command(arguments)
    return input_args_command

def make_save_shell_state_command():
    ## Save the shell state
    ## ```
    ## source $PASH_TOP/compiler/orchestration_runtime/save_shell_state.sh
    ## ```
    arguments = [string_to_argument("source"),
                 string_to_argument(config.SAVE_SHELL_STATE_EXECUTABLE)]
    input_args_command = make_command(arguments)
    return input_args_command

def restore_arguments_command():
    ## Restore the arguments to propagate internal changes, e.g., from `shift` outside.
    ## ```
    ## eval "set -- \"\${pash_input_args[@]}\""
    ## ```
    ##
    ## Alternative Solution: (TODO if we need extra performance -- avoiding eval) 
    ## Implement an AST node that accepts and returns a literal string
    ## bypassing unparsing. This would make this simpler and also more
    ## efficient (avoiding eval).
    ## However, it would require some work because we would need to implement
    ## support for this node in various places of PaSh and the unparser.
    ##      
    ##
    ## TODO: Maybe we need to only do this if there is a change.
    ## 
    set_arguments = [string_to_argument("eval"),
                     [['Q', string_to_argument('set -- ') +
                            [escaped_char('"')] + # The escaped quote
                            string_to_argument('\\${pash_input_args[@]}') +
                            [escaped_char('"')]]]]
    set_args_node = make_command(set_arguments)
    return set_args_node

def restore_exit_code_node():
    ## Restore the exit code (since now we have executed `set` last)
    ## ```
    ## ( exit "$pash_runtime_final_status")
    ## ```
    set_exit_status_command_arguments = [string_to_argument("exit"),
                                         [make_quoted_variable("pash_runtime_final_status")]]
    set_exit_status_command = make_command(set_exit_status_command_arguments)
    set_exit_status_node = make_kv('Subshell', [0, set_exit_status_command, []])
    return set_exit_status_node

## This function makes a command that calls the pash runtime
## together with the name of the file containing an IR. Then the
## pash runtime should read from this file and continue
## execution.
##
## TODO: At the moment this is written in python but it is in essense a simple shell script.
##       Is it possible to make it be a simple string instead of manually creating the AST?
##
## (MAYBE) TODO: The way I did it, is by calling the parser once, and seeing
## what it returns. Maybe it would make sense to call the parser on
## the fly to have a cleaner implementation here?
def make_call_to_pash_runtime(ir_filename, sequential_script_file_name,
                              disable_parallel_pipelines) -> AstNode:

    ## Disable parallel pipelines if we are in the last command of the script.
    ## ```
    ## pash_disable_parallel_pipelines=1
    ## ```
    if(disable_parallel_pipelines):
        assignments = [["pash_disable_parallel_pipelines",
                        string_to_argument("1")]]
    else:
        assignments = [["pash_disable_parallel_pipelines",
                        string_to_argument("0")]]
    assignments.append(["pash_sequential_script_file", 
                        string_to_argument(sequential_script_file_name)])
    assignments.append(["pash_input_ir_file", 
                        string_to_argument(ir_filename)])
    disable_parallel_pipelines_command = make_command([],
                                                      assignments=assignments)

    ## Call the runtime
    arguments = [string_to_argument("source"),
                 string_to_argument(config.RUNTIME_EXECUTABLE)]
    # ,
    #              string_to_argument(sequential_script_file_name),
    #              string_to_argument(ir_filename)]
    ## Pass all relevant argument to the planner
    ## TODO: Remove those
    # common_arguments_strings = config.pass_common_arguments(config.pash_args)
    # arguments += [string_to_argument(string) for string in common_arguments_strings]
    runtime_node = make_command(arguments)

    ## Create generic wrapper commands
    pre_runtime_nodes = make_pre_runtime_nodes()
    post_runtime_nodes = make_post_runtime_nodes()
    nodes = pre_runtime_nodes + [disable_parallel_pipelines_command, runtime_node] + post_runtime_nodes
    sequence = make_semi_sequence(nodes)
    return sequence

## TODO: Make that an actual call to the spec runtime
def make_call_to_spec_runtime(command_id: str) -> AstNode:
    ## Call the runtime
    arguments = [string_to_argument("source"),
                 string_to_argument(config.RUNTIME_EXECUTABLE),
                 string_to_argument(str(command_id))]
    ## Pass all relevant argument to the planner
    common_arguments_strings = config.pass_common_arguments(config.pash_args)
    arguments += [string_to_argument(string) for string in common_arguments_strings]
    runtime_node = make_command(arguments)

    ## Create generic wrapper commands
    pre_runtime_nodes = make_pre_runtime_nodes()
    post_runtime_nodes = make_post_runtime_nodes()
    nodes = pre_runtime_nodes + [runtime_node] + post_runtime_nodes
    sequence = make_semi_sequence(nodes)
    return sequence
