"""
AST to AST transformation

The preprocessing pass replaces all _candidate_ dataflow regions with
calls to PaSh's runtime to let it establish if they are actually dataflow
regions. The pass serializes all candidate dataflow regions:
- A list of ASTs if at the top level or
- an AST subtree if at a lower level

The PaSh runtime then deserializes them, compiles them (if safe) and optimizes them.
"""

from shasta.ast_node import AstNode

from shell_ast.ast_util import PreprocessedAST, UnparsedScript, unzip
from shell_ast.walk_preprocess import WalkPreprocess, PreprocessContext
from shell_ast.handlers.loop_tracking import for_node_with_loop_tracking
from shell_ast.transformation_options import AbstractTransformationState


# === Walker setup ===


def create_pash_walker() -> WalkPreprocess:
    """
    Create a preprocessing walker configured for PaSh.

    Returns:
        WalkPreprocess instance with PaSh-specific handlers
    """
    handlers = {
        "for": for_node_with_loop_tracking,
    }
    return WalkPreprocess(handlers=handlers)


# Module-level walker instance
_pash_walker = create_pash_walker()


def preprocess_node(
    ast_node: AstNode,
    trans_options: AbstractTransformationState,
    last_object: bool,
) -> PreprocessedAST:
    """
    Preprocesses an AstNode. Given an AstNode of any type, it will appropriately
    dispatch a preprocessor for the specific node type.

    This is the main entry point for preprocessing.

    Parameters:
        ast_node (AstNode): The AstNode to parse
        trans_options (AbstractTransformationState):
            A concrete transformation state instance corresponding to the output target
        last_object (bool): Flag for whether this is the last AstNode

    Returns:
        PreprocessedAst: the preprocessed version of the original AstNode
    """
    ctx = PreprocessContext(trans_options=trans_options, last_object=last_object)
    return _pash_walker.walk(ast_node, ctx)


# === Main transformation function ===


def replace_ast_regions(ast_objects, trans_options: AbstractTransformationState):
    """
    Replace candidate dataflow AST regions with calls to PaSh's runtime.
    """
    preprocessed_asts = []
    candidate_dataflow_region = []
    last_object = False
    for i, ast_object in enumerate(ast_objects):
        # If we are working on the last object we need to keep that in mind when replacing.
        # The last df-region should not be executed in parallel no matter what (to not lose its exit code.)
        if i == len(ast_objects) - 1:
            last_object = True

        ast, original_text, _linno_before, _linno_after = ast_object
        assert isinstance(ast, AstNode)

        # Preprocess ast by replacing subtrees with calls to runtime.
        # - If the whole AST needs to be replaced (e.g. if it is a pipeline)
        #   then the second output is true.
        # - If the next AST needs to be replaced too (e.g. if the current one is a background)
        #   then the third output is true
        preprocessed_ast_object = preprocess_node(
            ast, trans_options, last_object=last_object
        )
        # If the dataflow region is not maximal then it implies that the whole
        # AST should be replaced.
        assert (
            not preprocessed_ast_object.is_non_maximal()
            or preprocessed_ast_object.should_replace_whole_ast()
        )

        # If the whole AST needs to be replaced then it implies that
        # something will be replaced
        assert (
            not preprocessed_ast_object.should_replace_whole_ast()
            or preprocessed_ast_object.will_anything_be_replaced()
        )

        # If it isn't maximal then we just add it to the candidate
        if preprocessed_ast_object.is_non_maximal():
            candidate_dataflow_region.append(
                (preprocessed_ast_object.ast, original_text)
            )
        else:
            # If the current candidate dataflow region is non-empty
            # it means that the previous AST was in the background so
            # the current one has to be included in the process no matter what
            if len(candidate_dataflow_region) > 0:
                candidate_dataflow_region.append(
                    (preprocessed_ast_object.ast, original_text)
                )
                # Since the current one is maximal (or not wholy replaced)
                # we close the candidate.
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
                    # In this case, it is possible that no replacement happened,
                    # meaning that we can simply return the original parsed text as it was.
                    if (
                        preprocessed_ast_object.will_anything_be_replaced()
                        or original_text is None
                    ):
                        preprocessed_asts.append(preprocessed_ast_object.ast)
                    else:
                        preprocessed_asts.append(UnparsedScript(original_text))

    # Close the final dataflow region
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


def join_original_text_lines(shell_source_lines_or_none):
    """
    Join original unparsed shell source in a safe way,
    handling the case where some of the text is None (e.g., in case of stdin parsing).
    """
    if any([text_or_none is None for text_or_none in shell_source_lines_or_none]):
        return None
    else:
        return "\n".join(shell_source_lines_or_none)
