import copy

from shell_ast.ast_util import *
from shell_ast.transformation_options import AbstractTransformationState
from shasta.ast_node import AstNode


def preprocess_node(
    ast_node: AstNode,
    trans_options: AbstractTransformationState,
    last_object: bool,
) -> PreprocessedAST:
    """
    Preprocesses an AstNode. Given an AstNode of any type, it will appropriately
    dispatch a preprocessor for the specificy node type

    Parameters:
        ast_node (AstNode): The AstNode to parse
        trans_options (AbstractTransformationState):
            A concrete transformation state instance corresponding to the output target
        last_object (bool): Flag for whether this is the last AstNode

    Returns:
        PreprocessedAst: the preprocessed version of the original AstNode

    Note:
        For preprocess_node to dispatch the right function, the function being
        called must follow the convention "preprocess_node_<node_name>"
    """
    node_name = type(ast_node).NodeName.lower()
    preprocess_fn = globals().get(f"preprocess_node_{node_name}")
    if preprocess_fn is None:
        raise KeyError(f"Could not find appropriate preprocessor for {node_name}")
    return preprocess_fn(ast_node, trans_options, last_object)


## This preprocesses the AST node and also replaces it if it needs replacement .
## It is called by constructs that cannot be included in a dataflow region.
def preprocess_close_node(
    ast_node: AstNode,
    trans_options: AbstractTransformationState,
    last_object: bool = False,
):
    preprocessed_ast_object = preprocess_node(
        ast_node, trans_options, last_object=last_object
    )
    preprocessed_ast = preprocessed_ast_object.ast
    should_replace_whole_ast = preprocessed_ast_object.should_replace_whole_ast()
    if should_replace_whole_ast:
        final_ast = trans_options.replace_df_region(
            asts=[preprocessed_ast], disable_parallel_pipelines=last_object
        )
        something_replaced = True
    else:
        final_ast = preprocessed_ast
        something_replaced = preprocessed_ast_object.will_anything_be_replaced()
    return final_ast, something_replaced


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


def preprocess_node_pipe(
    ast_node: PipeNode,
    trans_options: AbstractTransformationState,
    last_object: bool = False,
):
    ## A pipeline is *always* a candidate dataflow region.
    ## Q: Is that true?

    ## TODO: Preprocess the internals of the pipe to allow
    ##       for mutually recursive calls to PaSh.
    ##
    ##       For example, if a command in the pipe has a command substitution
    ##       in one of its arguments then we would like to call our runtime
    ##       there instead of
    preprocessed_ast_object = PreprocessedAST(
        ast_node,
        replace_whole=True,
        non_maximal=ast_node.is_background,
        last_ast=last_object,
    )
    return preprocessed_ast_object


## TODO: Complete this
def preprocess_node_command(
    ast_node: CommandNode,
    trans_options: AbstractTransformationState,
    last_object: bool = False,
):
    ## TODO: Preprocess the internals of the pipe to allow
    ##       for mutually recursive calls to PaSh.
    ##
    ##       For example, if a command in the pipe has a command substitution
    ##       in one of its arguments then we would like to call our runtime
    ##       there instead of

    ## If there are no arguments, the command is just an
    ## assignment (Q: or just redirections?)
    if len(ast_node.arguments) == 0:
        preprocessed_ast_object = PreprocessedAST(
            ast_node,
            replace_whole=False,
            non_maximal=False,
            something_replaced=False,
            last_ast=last_object,
        )
        return preprocessed_ast_object

    ## This means we have a command. Commands are always candidate dataflow
    ## regions.
    preprocessed_ast_object = PreprocessedAST(
        ast_node, replace_whole=True, non_maximal=False, last_ast=last_object
    )
    return preprocessed_ast_object


# Background of (linno * t * redirection list)
## TODO: It might be possible to actually not close the inner node but rather apply the redirections on it
def preprocess_node_redir(
    ast_node: RedirNode,
    trans_options: AbstractTransformationState,
    last_object: bool = False,
):
    preprocessed_node, something_replaced = preprocess_close_node(
        ast_node.node, trans_options, last_object=last_object
    )
    ast_node.node = preprocessed_node
    preprocessed_ast_object = PreprocessedAST(
        ast_node,
        replace_whole=False,
        non_maximal=False,
        something_replaced=something_replaced,
        last_ast=last_object,
    )
    return preprocessed_ast_object


## TODO: Is that correct? Also, this should probably affect `semi`, `and`, and `or`
def preprocess_node_background(
    ast_node: BackgroundNode,
    trans_options: AbstractTransformationState,
    last_object: bool = False,
):
    ## A background node is *always* a candidate dataflow region.
    ## Q: Is that true?

    ## TODO: Preprocess the internals of the background to allow
    ##       for mutually recursive calls to PaSh.
    preprocessed_ast_object = PreprocessedAST(
        ast_node, replace_whole=True, non_maximal=True, last_ast=last_object
    )
    return preprocessed_ast_object


## TODO: We can actually preprocess the underlying node and then
##       return its characteristics above. However, we would need
##       to add a field in the IR that a node runs in a subshell
##       (which would have implications on how the backend outputs it).
##
##       e.g. a subshell node should also be output as a subshell in the backend.
## FIXME: This might not just be suboptimal, but also wrong.
def preprocess_node_subshell(
    ast_node: SubshellNode,
    trans_options: AbstractTransformationState,
    last_object: bool = False,
):
    preprocessed_body, something_replaced = preprocess_close_node(
        ast_node.body, trans_options, last_object=last_object
    )
    ast_node.body = preprocessed_body
    preprocessed_ast_object = PreprocessedAST(
        ast_node,
        replace_whole=False,
        non_maximal=False,
        something_replaced=something_replaced,
        last_ast=last_object,
    )
    return preprocessed_ast_object


## TODO: For all of the constructs below, think whether we are being too conservative


## TODO: This is not efficient at all since it calls the PaSh runtime everytime the loop is entered.
##       We have to find a way to improve that.
def preprocess_node_for(
    ast_node: ForNode,
    trans_options: AbstractTransformationState,
    last_object: bool = False,
):
    ## If we are in a loop, we push the loop identifier into the loop context
    loop_id = trans_options.enter_loop()
    preprocessed_body, something_replaced = preprocess_close_node(
        ast_node.body, trans_options, last_object=last_object
    )

    ## TODO: Then send this iteration identifier when talking to the spec scheduler
    ## TODO: After running checks put this behind a check to only run under speculation

    ## Create a new variable that tracks loop iterations
    var_name = loop_iter_var(loop_id)
    export_node = make_export_var_constant_string(var_name, "0")
    increment_node = make_increment_var(var_name)

    ## Also store the whole sequence of loop iters in a file
    all_loop_ids = trans_options.get_current_loop_context()

    ## export pash_loop_iters="$pash_loop_XXX_iter $pash_loop_YYY_iter ..."
    save_loop_iters_node = export_pash_loop_iters_for_current_context(all_loop_ids)

    ## Prepend the increment in the body
    ast_node.body = make_typed_semi_sequence(
        [
            to_ast_node(increment_node),
            to_ast_node(save_loop_iters_node),
            copy.deepcopy(preprocessed_body),
        ]
    )

    ## We pop the loop identifier from the loop context.
    ##
    ## KK 2023-04-27: Could this exit happen before the replacement leading to wrong
    ##     results? I think not because we use the _close_node preprocessing variant.
    ##     A similar issue might happen for while
    trans_options.exit_loop()

    ## reset the loop iters after we exit the loop
    out_of_loop_loop_ids = trans_options.get_current_loop_context()
    reset_loop_iters_node = export_pash_loop_iters_for_current_context(
        out_of_loop_loop_ids
    )

    ## Prepend the export in front of the loop
    # new_node = ast_node
    new_node = make_typed_semi_sequence(
        [to_ast_node(export_node), ast_node, to_ast_node(reset_loop_iters_node)]
    )
    # print(new_node)

    preprocessed_ast_object = PreprocessedAST(
        new_node,
        replace_whole=False,
        non_maximal=False,
        something_replaced=something_replaced,
        last_ast=last_object,
    )

    return preprocessed_ast_object


def preprocess_node_while(
    ast_node: WhileNode,
    trans_options: AbstractTransformationState,
    last_object: bool = False,
):
    ## If we are in a loop, we push the loop identifier into the loop context
    trans_options.enter_loop()

    preprocessed_test, sth_replaced_test = preprocess_close_node(
        ast_node.test, trans_options, last_object=last_object
    )
    preprocessed_body, sth_replaced_body = preprocess_close_node(
        ast_node.body, trans_options, last_object=last_object
    )
    ast_node.test = preprocessed_test
    ast_node.body = preprocessed_body
    something_replaced = sth_replaced_test or sth_replaced_body
    preprocessed_ast_object = PreprocessedAST(
        ast_node,
        replace_whole=False,
        non_maximal=False,
        something_replaced=something_replaced,
        last_ast=last_object,
    )

    ## We pop the loop identifier from the loop context.
    trans_options.exit_loop()
    return preprocessed_ast_object


## This is the same as the one for `For`
def preprocess_node_defun(
    ast_node: DefunNode,
    trans_options: AbstractTransformationState,
    last_object: bool = False,
):
    ## TODO: For now we don't want to compile function bodies
    # preprocessed_body = preprocess_close_node(ast_node.body)
    # ast_node.body = preprocessed_body
    preprocessed_ast_object = PreprocessedAST(
        ast_node,
        replace_whole=False,
        non_maximal=False,
        something_replaced=False,
        last_ast=last_object,
    )
    return preprocessed_ast_object


## TODO: If the preprocessed is not maximal we actually need to combine it with the one on the right.
def preprocess_node_semi(
    ast_node: SemiNode,
    trans_options: AbstractTransformationState,
    last_object: bool = False,
):
    # preprocessed_left, should_replace_whole_ast, is_non_maximal = preprocess_node(ast_node.left, irFileGen, config)
    ##
    ## TODO: Is it valid that only the right one is considered the last command?
    preprocessed_left, sth_replaced_left = preprocess_close_node(
        ast_node.left_operand, trans_options, last_object
    )
    preprocessed_right, sth_replaced_right = preprocess_close_node(
        ast_node.right_operand, trans_options, last_object=last_object
    )
    ast_node.left_operand = preprocessed_left
    ast_node.right_operand = preprocessed_right
    sth_replaced = sth_replaced_left or sth_replaced_right
    preprocessed_ast_object = PreprocessedAST(
        ast_node,
        replace_whole=False,
        non_maximal=False,
        something_replaced=sth_replaced,
        last_ast=last_object,
    )
    return preprocessed_ast_object


## TODO: Make sure that what is inside an `&&`, `||`, `!` (and others) does not run in parallel_pipelines
##       since we need its exit code.
def preprocess_node_and(
    ast_node: AndNode,
    trans_options: AbstractTransformationState,
    last_object: bool = False,
):
    # preprocessed_left, should_replace_whole_ast, is_non_maximal = preprocess_node(ast_node.left, irFileGen, config)
    preprocessed_left, sth_replaced_left = preprocess_close_node(
        ast_node.left_operand, trans_options, last_object=last_object
    )
    preprocessed_right, sth_replaced_right = preprocess_close_node(
        ast_node.right_operand, trans_options, last_object=last_object
    )
    ast_node.left_operand = preprocessed_left
    ast_node.right_operand = preprocessed_right
    sth_replaced = sth_replaced_left or sth_replaced_right
    preprocessed_ast_object = PreprocessedAST(
        ast_node,
        replace_whole=False,
        non_maximal=False,
        something_replaced=sth_replaced,
        last_ast=last_object,
    )
    return preprocessed_ast_object


def preprocess_node_or(
    ast_node: OrNode,
    trans_options: AbstractTransformationState,
    last_object: bool = False,
):
    # preprocessed_left, should_replace_whole_ast, is_non_maximal = preprocess_node(ast_node.left, irFileGen, config)
    preprocessed_left, sth_replaced_left = preprocess_close_node(
        ast_node.left_operand, trans_options, last_object=last_object
    )
    preprocessed_right, sth_replaced_right = preprocess_close_node(
        ast_node.right_operand, trans_options, last_object=last_object
    )
    ast_node.left_operand = preprocessed_left
    ast_node.right_operand = preprocessed_right
    sth_replaced = sth_replaced_left or sth_replaced_right
    preprocessed_ast_object = PreprocessedAST(
        ast_node,
        replace_whole=False,
        non_maximal=False,
        something_replaced=sth_replaced,
        last_ast=last_object,
    )
    return preprocessed_ast_object


def preprocess_node_not(
    ast_node: NotNode,
    trans_options: AbstractTransformationState,
    last_object: bool = False,
):
    # preprocessed_left, should_replace_whole_ast, is_non_maximal = preprocess_node(ast_node.left)
    preprocessed_body, sth_replaced = preprocess_close_node(
        ast_node.body, trans_options, last_object=last_object
    )
    ast_node.body = preprocessed_body
    preprocessed_ast_object = PreprocessedAST(
        ast_node,
        replace_whole=False,
        non_maximal=False,
        something_replaced=sth_replaced,
        last_ast=last_object,
    )
    return preprocessed_ast_object


def preprocess_node_if(
    ast_node: IfNode,
    trans_options: AbstractTransformationState,
    last_object: bool = False,
):
    # preprocessed_left, should_replace_whole_ast, is_non_maximal = preprocess_node(ast_node.left, irFileGen, config)
    preprocessed_cond, sth_replaced_cond = preprocess_close_node(
        ast_node.cond, trans_options, last_object=last_object
    )
    preprocessed_then, sth_replaced_then = preprocess_close_node(
        ast_node.then_b, trans_options, last_object=last_object
    )
    preprocessed_else, sth_replaced_else = preprocess_close_node(
        ast_node.else_b, trans_options, last_object=last_object
    )
    ast_node.cond = preprocessed_cond
    ast_node.then_b = preprocessed_then
    ast_node.else_b = preprocessed_else
    sth_replaced = sth_replaced_cond or sth_replaced_then or sth_replaced_else
    preprocessed_ast_object = PreprocessedAST(
        ast_node,
        replace_whole=False,
        non_maximal=False,
        something_replaced=sth_replaced,
        last_ast=last_object,
    )
    return preprocessed_ast_object


def preprocess_case(
    case, trans_options: AbstractTransformationState, last_object: bool
):
    preprocessed_body, sth_replaced = preprocess_close_node(
        case["cbody"], trans_options, last_object=last_object
    )
    case["cbody"] = preprocessed_body
    return case, sth_replaced


def preprocess_node_case(
    ast_node: CaseNode,
    trans_options: AbstractTransformationState,
    last_object: bool = False,
):
    preprocessed_cases_replaced = [
        preprocess_case(case, trans_options, last_object=last_object)
        for case in ast_node.cases
    ]
    preprocessed_cases, sth_replaced_cases = list(zip(*preprocessed_cases_replaced))
    ast_node.cases = preprocessed_cases
    preprocessed_ast_object = PreprocessedAST(
        ast_node,
        replace_whole=False,
        non_maximal=False,
        something_replaced=any(sth_replaced_cases),
        last_ast=last_object,
    )
    return preprocessed_ast_object
