"""
AST to AST transformation


The preprocessing pass replaces all _candidate_ dataflow regions with
calls to PaSh's runtime to let it establish if they are actually dataflow
regions. The pass serializes all candidate dataflow regions:
- A list of ASTs if at the top level or
- an AST subtree if at a lower level

The PaSh runtime then deserializes the(m, compiles them (if safe) and optimizes them.
"""

from env_var_names import *
from shell_ast.ast_util import *
from shell_ast.preprocess_ast_cases import preprocess_node
from shell_ast.transformation_options import AbstractTransformationState


def replace_ast_regions(ast_objects, trans_options: AbstractTransformationState):
    """
    Replace candidate dataflow AST regions with calls to PaSh's runtime.
    """
    preprocessed_asts = []
    candidate_dataflow_region = []
    last_object = False
    for i, ast_object in enumerate(ast_objects):
        # log("Preprocessing AST {}".format(i))
        # log(ast_object)
        ## If we are working on the last object we need to keep that in mind when replacing.
        ##
        ## The last df-region should not be executed in parallel no matter what (to not lose its exit code.)
        if i == len(ast_objects) - 1:
            # log("Last object")
            last_object = True

        ast, original_text, _linno_before, _linno_after = ast_object
        assert isinstance(ast, AstNode)

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
        preprocessed_ast_object = preprocess_node(
            ast, trans_options, last_object=last_object
        )
        ## If the dataflow region is not maximal then it implies that the whole
        ## AST should be replaced.
        assert (
            not preprocessed_ast_object.is_non_maximal()
            or preprocessed_ast_object.should_replace_whole_ast()
        )

        ## If the whole AST needs to be replaced then it implies that
        ## something will be replaced
        assert (
            not preprocessed_ast_object.should_replace_whole_ast()
            or preprocessed_ast_object.will_anything_be_replaced()
        )

        ## If it isn't maximal then we just add it to the candidate
        if preprocessed_ast_object.is_non_maximal():
            candidate_dataflow_region.append(
                (preprocessed_ast_object.ast, original_text)
            )
        else:
            ## If the current candidate dataflow region is non-empty
            ## it means that the previous AST was in the background so
            ## the current one has to be included in the process no matter what
            if len(candidate_dataflow_region) > 0:
                candidate_dataflow_region.append(
                    (preprocessed_ast_object.ast, original_text)
                )
                ## Since the current one is maximal (or not wholy replaced)
                ## we close the candidate.
                dataflow_region_asts, dataflow_region_lines = unzip(
                    candidate_dataflow_region
                )
                dataflow_region_text = join_original_text_lines(dataflow_region_lines)
                replaced_ast = trans_options.replace_df_region(
                    dataflow_region_asts,
                    ast_text=dataflow_region_text,
                    disable_parallel_pipelines=last_object,
                )
                candidate_dataflow_region = []
                preprocessed_asts.append(replaced_ast)
            else:
                if preprocessed_ast_object.should_replace_whole_ast():
                    replaced_ast = trans_options.replace_df_region(
                        [preprocessed_ast_object.ast],
                        ast_text=original_text,
                        disable_parallel_pipelines=last_object,
                    )
                    preprocessed_asts.append(replaced_ast)
                else:
                    ## In this case, it is possible that no replacement happened,
                    ## meaning that we can simply return the original parsed text as it was.
                    if (
                        preprocessed_ast_object.will_anything_be_replaced()
                        or original_text is None
                    ):
                        preprocessed_asts.append(preprocessed_ast_object.ast)
                    else:
                        preprocessed_asts.append(UnparsedScript(original_text))

    ## Close the final dataflow region
    if len(candidate_dataflow_region) > 0:
        dataflow_region_asts, dataflow_region_lines = unzip(candidate_dataflow_region)
        dataflow_region_text = join_original_text_lines(dataflow_region_lines)
        replaced_ast = trans_options.replace_df_region(
            dataflow_region_asts,
            ast_text=dataflow_region_text,
            disable_parallel_pipelines=True,
        )
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
