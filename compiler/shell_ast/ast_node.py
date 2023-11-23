"""
Wrapper classes for shasta's Commands.
The motivation for creating wrappers around shasta's Commands are that
we want to add PaSh specific methods to these Commands rather than doing an
unsave and expensive typecheck
"""

from abc import ABC, abstractmethod
from copy import deepcopy
from typing import Any, Protocol, Tuple, TypedDict, TYPE_CHECKING
from pash_annotations.annotation_generation.datatypes.InputOutputInfo import InputOutputInfo
from pash_annotations.datatypes.CommandInvocationInitial import CommandInvocationInitial

import shasta.ast_node
from shasta.ast_node import make_typed_semi_sequence
from shasta.json_to_ast import to_ast_node

from shell_ast.ast_util import (
    PreprocessedAST,
    loop_iter_var,
    make_export_var_constant_string,
    make_increment_var,
    export_pash_loop_iters_for_current_context,
)
from ir import IR, compile_command_to_DFG
from util import log

if TYPE_CHECKING:
    from transformation_options import AbstractTransformationState as TransState


class PashNode(Protocol):
    """
    A protocol class representing the PaSh wrapper classes for shasta's Command.
    It checks that a class has all methods that a shasta Command would have, plus
    the PashMixin methods. It is somewhat analogous to an interface.

    Because the concrete PashNodes inherit from a concrete subclass of Command rather
    than Command itself, a Protocol is required to adequately represent it.
    As outlined in the Python docs, this is mostly for the sake of accurate
    type-hinting, and does not have any runtime implications
    """

    def json(self):
        raise NotImplementedError

    def pretty(self) -> str:
        raise NotImplementedError

    def preprocess(
        self, trans_state: "TransState", last_object: bool
    ) -> PreprocessedAST:
        raise NotImplementedError

    def _preprocess_close_node(self, trans_state: "TransState", last_object: bool):
        raise NotImplementedError

    def compile(self, file_id_gen, config) -> IR:
        raise NotImplementedError


class PashMixin(ABC):
    """
    Mixin that extends functionality of shasta's Commands
    """

    def _preprocess_close_node(
        self: PashNode, trans_state: "TransState", last_object: bool = False
    ) -> Tuple[PashNode, bool]:
        preprocessed_ast_object = self.preprocess(trans_state, last_object=last_object)
        preprocessed_ast = preprocessed_ast_object.ast
        should_replace_whole_ast = preprocessed_ast_object.should_replace_whole_ast()

        if should_replace_whole_ast:
            final_ast = trans_state.replace_df_region(
                asts=[preprocessed_ast], disable_parallel_pipelines=last_object
            )
            something_replaced = True
        else:
            final_ast = preprocessed_ast
            something_replaced = preprocessed_ast_object.will_anything_be_replaced()
        return final_ast, something_replaced

    def _compile_arg_char(self, arg, file_id_gen, config):
        """
        Compile an arg char recursively if it contains quotes or command substitution.

        Currently being extended to also expand any arguments taht are safe to expand.
        """
        ## Compile the arg char
        if isinstance(arg, shasta.ast_node.CArgChar) or isinstance(
            arg, shasta.ast_node.EArgChar
        ):
            # Single character or escape
            return arg
        elif isinstance(arg, shasta.ast_node.BArgChar):
            # TODO: I probably have to redirect the input of the compiled
            #       node (IR) to be closed, and the output to be
            #       redirected to some file that we will use to write to
            #       the command argument to complete the command
            #       substitution.
            arg.node = arg.node.compile(arg.node, file_id_gen, config)
            return arg
        elif isinstance(arg, shasta.ast_node.QArgChar):
            arg.arg = self._compile_command_argument(arg.arg, file_id_gen, config)
            return arg
        else:
            log(f"Unknown arg_char: {arg}")
            ## TODO: Complete this
            return arg

    def _compile_command_to_dfg(file_id_gen, command_name, options, redirections):


    def _compile_command_argument(self, arg, file_id_gen, config):
        return [self._compile_arg_char(char, file_id_gen, config) for char in arg]

    def _compile_assignments(self, assignments, file_id_gen, config):
        return [
            [
                assignment[0],
                self._compile_command_argument(assignment[1], file_id_gen, config),
            ]
            for assignment in assignments
        ]

    def _compile_redirection(self, redir, file_id_gen, config):
        file_arg = self._compile_command_argument(redir.arg, file_id_gen, config)
        redir.arg = file_arg
        return redir


class PipeNode(shasta.ast_node.PipeNode, PashMixin):
    items: list[PashNode]

    def __init__(self, is_background, items):
        super().__init__(is_background, items)
        self.items = [pash_node_from(n) for n in items]

    def preprocess(self, trans_state: "TransState", last_object: bool = False):
        preprocessed_ast_object = PreprocessedAST(
            self,
            replace_whole=True,
            non_maximal=self.is_background,
            last_ast=last_object,
        )
        return preprocessed_ast_object

    def compile(self, file_id_gen, config):
        compiled_ir = self._combine_pipe(
            [item.compile(file_id_gen, config) for item in self.items]
        )
        compiled_ir.set_ast(self.json())
        compiled_ir.set_background(self.is_background)
        return compiled_ir

    def _combine_pipe(self, items):
        combined_nodes = items[0]
        for item in items[1:]:
            combined_nodes.pipe_append(item)

        return combined_nodes


class CommandNode(shasta.ast_node.CommandNode, PashMixin):
    def __init__(self, line_number, assignments, arguments, redir_list):
        super().__init__(line_number, assignments, arguments, redir_list)

    def preprocess(self, trans_state: "TransState", last_object: bool = False):
        if len(self.arguments) == 0:
            return PreprocessedAST(
                self,
                replace_whole=False,
                non_maximal=False,
                something_replaced=False,
                last_ast=last_object,
            )
        else:
            return PreprocessedAST(
                self, replace_whole=True, non_maximal=False, last_ast=last_object
            )

    def compile(self, file_id_gen, config):
        compiled_assignments = [
            self._compile_assignments(self.assignments, file_id_gen, config)
        ]
        compiled_redirs = [
            self._compile_redirection(redir, file_id_gen, config)
            for redir in self.redir_list
        ]

        assert len(ast_node.arguments) > 0

        args = self.arguments
        command_name = args[0]
        options = [
            self._compile_command_argument(arg, file_id_gen, config) for arg in args[1:]
        ]

        try:
            ir = compile_command_to_DFG(file_id_gen, command_name, options, compiled_redirs)
            compiled_ast = ir
        except ValueError as err:
            log("Command not compiled to DFG:", err)
            ## TODO: Maybe we want to fail here instead of waiting for later?
            ##       Is there any case where a non-compiled command is fine?
            # log(traceback.format_exc())
            compiled_arguments = compile_command_arguments(arguments, fileIdGen, config)
            compiled_ast = make_kv(
                type(ast_node).NodeName,
                [
                    ast_node.line_number,
                    compiled_assignments,
                    compiled_arguments,
                    compiled_redirections,
                ],
            )

        return compiled_ast


class SubshellNode(shasta.ast_node.SubshellNode, PashMixin):
    body: PashNode

    def __init__(self, line_number, body, redir_list):
        super().__init__(line_number, body, redir_list)
        self.body = pash_node_from(body)

    def preprocess(self, trans_state: "TransState", last_object: bool = False):
        preprocessed_body, something_replaced = self._preprocess_close_node(
            trans_state, last_object=last_object
        )
        self.body = preprocessed_body
        preprocessed_ast_object = PreprocessedAST(
            self,
            replace_whole=False,
            non_maximal=False,
            something_replaced=something_replaced,
            last_ast=last_object,
        )
        return preprocessed_ast_object


class AndNode(shasta.ast_node.AndNode, PashMixin):
    left_operand: PashNode
    right_operand: PashNode

    def __init__(self, left_operand, right_operand):
        super().__init__(left_operand, right_operand)
        self.left_operand = pash_node_from(left_operand)
        self.right_operand = pash_node_from(right_operand)

    def preprocess(self, trans_state: "TransState", last_object: bool = False):
        preprocessed_l, sth_replaced_l = self.left_operand(
            trans_state, last_object=last_object
        )
        preprocessed_r, sth_replaced_r = self.right_operand(
            trans_state, last_object=last_object
        )
        self.left_operand = preprocessed_l
        self.right_operand = preprocessed_r
        sth_replaced = sth_replaced_l or sth_replaced_r
        preprocessed_ast_object = PreprocessedAST(
            self,
            replace_whole=False,
            non_maximal=False,
            something_replaced=sth_replaced,
            last_ast=last_object,
        )
        return preprocessed_ast_object


class OrNode(shasta.ast_node.OrNode, PashMixin):
    left_operand: PashNode
    right_operand: PashNode

    def __init__(self, left_operand, right_operand):
        super().__init__(left_operand, right_operand)
        self.left_operand = pash_node_from(left_operand)
        self.right_operand = pash_node_from(right_operand)

    def preprocess(
        self, trans_state: "TransState", last_object: bool = False
    ) -> PreprocessedAST:
        preprocessed_l, sth_replaced_l = self.left_operand(
            trans_state, last_object=last_object
        )
        preprocessed_r, sth_replaced_r = self.right_operand(
            trans_state, last_object=last_object
        )
        self.left_operand = preprocessed_l
        self.right_operand = preprocessed_r
        sth_replaced = sth_replaced_l or sth_replaced_r
        preprocessed_ast_object = PreprocessedAST(
            self,
            replace_whole=False,
            non_maximal=False,
            something_replaced=sth_replaced,
            last_ast=last_object,
        )
        return preprocessed_ast_object


# TODO: If the preprocessed is not maximal we actually need to combine it with the one on the right.
class SemiNode(shasta.ast_node.SemiNode, PashMixin):
    left_operand: PashNode
    right_operand: PashNode

    def __init__(self, left_operand, right_operand):
        super().__init__(left_operand, right_operand)
        self.left_operand = pash_node_from(left_operand)
        self.right_operand = pash_node_from(right_operand)

    def preprocess(
        self, trans_state: "TransState", last_object: bool
    ) -> PreprocessedAST:
        # TODO: Is it valid that only the right one is considered the last command?
        preprocessed_l, sth_replaced_l = self.left_operand(trans_state, last_object)
        preprocessed_r, sth_replaced_r = self.right_operand(trans_state, last_object)
        self.left_operand = preprocessed_l
        self.right_operand = preprocessed_r
        sth_replaced = sth_replaced_l or sth_replaced_r
        return PreprocessedAST(
            self,
            replace_whole=False,
            non_maximal=False,
            something_replaced=sth_replaced,
            last_ast=last_object,
        )


class NotNode(shasta.ast_node.NotNode, PashMixin):
    body: PashNode

    def __init__(self, body):
        super().__init__(body)
        self.body = pash_node_from(body)

    def preprocess(
        self, trans_state: "TransState", last_object: bool = False
    ) -> PreprocessedAST:
        preprocessed_body, sth_replaced = self.body(
            trans_state, last_object=last_object
        )
        self.body = preprocessed_body
        return PreprocessedAST(
            self,
            replace_whole=False,
            non_maximal=False,
            something_replaced=sth_replaced,
            last_ast=last_object,
        )


class RedirNode(shasta.ast_node.RedirNode, PashMixin):
    line_number: int
    node: PashNode
    redir_list: list

    def __init__(self, line_number, node, redir_list):
        super().__init__(line_number, node, redir_list)
        self.node = pash_node_from(node)

    def preprocess(
        self, trans_state: "TransState", last_object: bool
    ) -> PreprocessedAST:
        preprocessed_node, sth_replaced = self.node(trans_state, last_object)
        self.node = preprocessed_node
        return PreprocessedAST(
            self,
            replace_whole=False,
            non_maximal=False,
            something_replaced=sth_replaced,
            last_ast=last_object,
        )


class BackgroundNode(shasta.ast_node.BackgroundNode, PashMixin):
    node: PashNode

    def __init__(self, line_number, node, redir_list):
        super().__init__(line_number, node, redir_list)
        self.node = pash_node_from(node)

    def preprocess(
        self, trans_state: "TransState", last_object: bool
    ) -> PreprocessedAST:
        return PreprocessedAST(
            self, replace_whole=True, non_maximal=True, last_ast=last_object
        )


class DefunNode(shasta.ast_node.DefunNode, PashMixin):
    body: PashNode

    def __init__(self, line_number, name, body):
        super().__init__(line_number, name, body)
        self.body = pash_node_from(body)

    def preprocess(
        self, trans_state: "TransState", last_object: bool
    ) -> PreprocessedAST:
        ## TODO: For now we don't want to compile function bodies
        # preprocessed_body = preprocess_close_node(ast_node.body)
        # ast_node.body = preprocessed_body
        return PreprocessedAST(
            self,
            replace_whole=False,
            non_maximal=False,
            something_replaced=False,
            last_ast=last_object,
        )


class ForNode(shasta.ast_node.ForNode, PashMixin):
    body: PashNode

    def __init__(self, line_number, argument, body, variable):
        super().__init__(line_number, argument, body, variable)
        self.body = pash_node_from(body)

    def preprocess(
        self, trans_state: "TransState", last_object: bool
    ) -> PreprocessedAST:
        # If we are in a loop, we push the loop identifier into the loop context
        loop_id = trans_state.enter_loop()
        preprocessed_body, sth_replaced = self.body._preprocess_close_node(
            trans_state, last_object
        )
        # TODO: Then send this iteration identifier when talking to the spec scheduler
        # TODO: After running checks put this behind a check to only run under speculation

        # Create a new variable that tracks loop iterations
        var_name = loop_iter_var(loop_id)
        export_node = make_export_var_constant_string(var_name, "0")
        increment_node = make_increment_var(var_name)

        # Also store the whole sequence of loop iters in a file
        all_loop_ids = trans_state.get_current_loop_context()

        # export pash_loop_iters="$pash_loop_XXX_iter $pash_loop_YYY_iter ..."
        save_loop_iters_node = export_pash_loop_iters_for_current_context(all_loop_ids)

        self.body = pash_node_from(
            make_typed_semi_sequence(
                [
                    to_ast_node(increment_node),
                    to_ast_node(save_loop_iters_node),
                    deepcopy(preprocessed_body),
                ]
            )
        )

        # We pop the loop identifier from the loop context.
        #
        # KK 2023-04-27: Could this exit happen before the replacement leading to wrong
        #     results? I think not because we use the _close_node preprocessing variant.
        #     A similar issue might happen for while
        trans_state.exit_loop()

        # reset the loop iters after we exit the loop
        out_of_loop_loop_ids = trans_state.get_current_loop_context()
        reset_loop_iters_node = export_pash_loop_iters_for_current_context(
            out_of_loop_loop_ids
        )

        # Prepend the export in front of the loop
        # new_node = ast_node
        new_node = pash_node_from(
            make_typed_semi_sequence(
                [to_ast_node(export_node), self, to_ast_node(reset_loop_iters_node)]
            )
        )
        # print(new_node)

        preprocessed_ast_object = PreprocessedAST(
            new_node,
            replace_whole=False,
            non_maximal=False,
            something_replaced=sth_replaced,
            last_ast=last_object,
        )

        return preprocessed_ast_object


class WhileNode(shasta.ast_node.WhileNode, PashMixin):
    test: PashNode
    body: PashNode

    def __init__(self, test, body):
        super().__init__(test, body)
        self.test = pash_node_from(test)
        self.body = pash_node_from(body)

    def preprocess(
        self, trans_state: "TransState", last_object: bool
    ) -> PreprocessedAST:
        trans_state.enter_loop()

        preprocessed_test, sth_replaced_test = self.test(trans_state, last_object)
        preprocessed_body, sth_replaced_body = self.body(trans_state, last_object)
        self.test = preprocessed_test
        self.body = preprocessed_body
        sth_replaced = sth_replaced_test or sth_replaced_body
        preprocessed_ast_object = PreprocessedAST(
            self,
            replace_whole=False,
            non_maximal=False,
            something_replaced=sth_replaced,
            last_ast=last_object,
        )
        trans_state.exit_loop()
        return preprocessed_ast_object


class IfNode(shasta.ast_node.IfNode, PashMixin):
    cond: PashNode
    then_b: PashNode
    else_b: PashNode

    def __init__(self, cond, then_b, else_b):
        super().__init__(cond, then_b, else_b)
        self.cond = pash_node_from(cond)
        self.then_b = pash_node_from(then_b)
        self.else_b = pash_node_from(else_b)

    def preprocess(
        self, trans_state: "TransState", last_object: bool
    ) -> PreprocessedAST:
        preprocess_cond, sth_replaced_cond = self.cond(trans_state, last_object)
        preprocess_then, sth_replaced_then = self.then_b(trans_state, last_object)
        preprocess_else, sth_replaced_else = self.else_b(trans_state, last_object)
        self.cond = preprocess_cond
        self.then_b = preprocess_then
        self.else_b = preprocess_else
        sth_replaced = sth_replaced_cond or sth_replaced_then or sth_replaced_else
        return PreprocessedAST(
            self,
            replace_whole=False,
            non_maximal=False,
            something_replaced=sth_replaced,
            last_ast=last_object,
        )


class CaseNode(shasta.ast_node.CaseNode, PashMixin):
    class CaseContent(TypedDict):
        cpattern: Any
        cbody: PashNode

    cases: list[CaseContent]

    def __init__(self, line_number, argument, cases):
        super().__init__(line_number, argument, cases)
        self.cases = []
        for case in cases:
            case["cbody"] = pash_node_from(case)
            self.cases.append(case)

    def preprocess(
        self, trans_state: "TransState", last_object: bool
    ) -> PreprocessedAST:
        preprocessed_cases_replaced = [
            self._preprocess_case(case, trans_state, last_object) for case in self.cases
        ]
        preprocessed_cases, sth_replaced_cases = list(zip(*preprocessed_cases_replaced))
        self.cases = preprocessed_cases
        return PreprocessedAST(
            self,
            replace_whole=False,
            non_maximal=False,
            something_replaced=any(sth_replaced_cases),
            last_ast=last_object,
        )

    def _preprocess_case(self, case, trans_state, last_object):
        preprocessed_body, sth_replaced = case["cbody"](trans_state, last_object)
        case["cbody"] = preprocessed_body
        return case, sth_replaced


def pash_node_from(shasta_node: shasta.ast_node.AstNode) -> PashNode:
    pash_node_class = globals()[f"{shasta_node.NodeName}Node"]
    args = vars(shasta_node)
    return pash_node_class(**args)
