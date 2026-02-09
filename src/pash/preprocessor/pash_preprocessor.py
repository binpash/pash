"""
PaSh Preprocessor

This module handles preprocessing of shell scripts by:
1. Parsing the shell script to ASTs
2. Replacing candidate dataflow regions with calls to PaSh runtime
3. Unparsing the transformed ASTs back to shell syntax
"""

import sys
import os
import argparse
import logging
from datetime import datetime

from shasta.ast_node import AstNode

from ast_util import PreprocessedAST, UnparsedScript, unzip
from transformation_options import TransformationState
from walk_preprocess import WalkPreprocess, PreprocessContext
from parse import parse_shell_to_asts, from_ast_objects_to_shell
from util import log, logging_prefix, print_time_delta


LOGGING_PREFIX = "PaSh Preprocessor: "

## Increase the recursion limit (it seems that the parser/unparser needs it for bigger scripts)
sys.setrecursionlimit(10000)
## Note: The preprocessor is very slow for very large recursive scripts


# Module-level walker instance
_pash_walker = WalkPreprocess()


def config_from_args(pash_args):
    """Configure logging based on command-line arguments."""
    if pash_args.log_file == "":
        logging.basicConfig(format="%(message)s")
    else:
        logging.basicConfig(
            format="%(message)s",
            filename=f"{os.path.abspath(pash_args.log_file)}",
            filemode="w",
        )

    if pash_args.debug == 0:
        logging.getLogger().setLevel(logging.ERROR)
    elif pash_args.debug == 1:
        logging.getLogger().setLevel(logging.WARNING)
    elif pash_args.debug == 2:
        logging.getLogger().setLevel(logging.INFO)
    elif pash_args.debug >= 3:
        logging.getLogger().setLevel(logging.DEBUG)


class Parser(argparse.ArgumentParser):
    """Command-line argument parser for the preprocessor."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.add_argument(
            "-d",
            "--debug",
            type=int,
            help="configure debug level; defaults to 0",
            default=1,
        )
        self.add_argument(
            "--log_file",
            help="configure where to write the log; defaults to stderr.",
            default="",
        )
        self.add_argument(
            "input",
            help="the script file to be preprocessed",
        )
        self.add_argument(
            "--output",
            help="path where the preprocessed script will be saved",
            required=True,
        )
        self.add_argument(
            "--bash",
            help="(experimental) interpret the input as a bash script file",
            action="store_true",
        )


@logging_prefix(LOGGING_PREFIX)
def main():
    """Main entry point for the preprocessor."""
    args = parse_args()
    input_script_path = args.input

    preprocessed_shell_script = preprocess(input_script_path, args)

    # Write the preprocessed script to the output file
    fname = args.output
    log("Preprocessed script stored in:", fname)
    with open(fname, "wb") as new_shell_file:
        preprocessed_shell_script = preprocessed_shell_script.encode(
            "utf-8", errors="replace"
        )
        new_shell_file.write(preprocessed_shell_script)

    log("-" * 40)  # log end marker


def preprocess(input_script_path, args):
    """
    Preprocess a shell script.

    This function:
    1. Parses the shell script to ASTs
    2. Preprocesses ASTs by replacing candidate dataflow regions
    3. Unparses the ASTs back to shell syntax

    Args:
        input_script_path: Path to the input shell script
        args: Parsed command-line arguments

    Returns:
        The preprocessed shell script as a string
    """
    # 1. Parse shell to AST
    preprocessing_parsing_start_time = datetime.now()
    ast_objects = parse_shell_to_asts(input_script_path, bash_mode=args.bash)
    preprocessing_parsing_end_time = datetime.now()
    print_time_delta(
        "Preprocessing -- Parsing",
        preprocessing_parsing_start_time,
        preprocessing_parsing_end_time,
    )

    # 2. Preprocess ASTs by replacing candidates with calls to PaSh runtime
    preprocessing_pash_start_time = datetime.now()
    preprocessed_asts = preprocess_asts(ast_objects, args)
    preprocessing_pash_end_time = datetime.now()
    print_time_delta(
        "Preprocessing -- PaSh",
        preprocessing_pash_start_time,
        preprocessing_pash_end_time,
    )

    # 3. Unparse the ASTs back to shell syntax
    preprocessing_unparsing_start_time = datetime.now()
    preprocessed_shell_script = from_ast_objects_to_shell(preprocessed_asts)
    preprocessing_unparsing_end_time = datetime.now()
    print_time_delta(
        "Preprocessing -- Unparsing",
        preprocessing_unparsing_start_time,
        preprocessing_unparsing_end_time,
    )

    return preprocessed_shell_script


def preprocess_asts(ast_objects, args):
    """
    Preprocess AST objects by replacing candidate dataflow regions.

    Args:
        ast_objects: List of parsed AST objects
        args: Parsed command-line arguments

    Returns:
        List of preprocessed AST objects
    """
    trans_options = TransformationState()
    return replace_ast_regions(ast_objects, trans_options)


# === AST region replacement ===


def preprocess_node(
    ast_node: AstNode,
    trans_options: TransformationState,
    last_object: bool,
) -> PreprocessedAST:
    """
    Preprocesses an AstNode. Given an AstNode of any type, it will appropriately
    dispatch a preprocessor for the specific node type.

    Parameters:
        ast_node (AstNode): The AstNode to parse
        trans_options (TransformationState):
            A concrete transformation state instance corresponding to the output target
        last_object (bool): Flag for whether this is the last AstNode

    Returns:
        PreprocessedAst: the preprocessed version of the original AstNode
    """
    ctx = PreprocessContext(trans_options=trans_options, last_object=last_object)
    return _pash_walker.walk(ast_node, ctx)


def replace_ast_regions(ast_objects, trans_options: TransformationState):
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


def parse_args():
    """Parse command-line arguments."""
    prog_name = sys.argv[0]
    if "PASH_FROM_SH" in os.environ:
        prog_name = os.environ["PASH_FROM_SH"]
    parser = Parser(prog_name)
    args = parser.parse_args()
    config_from_args(args)

    # Log all arguments
    log("Arguments:")
    for arg_name, arg_val in vars(args).items():
        log(arg_name, arg_val)
    log("-" * 40)

    return args


if __name__ == "__main__":
    main()
