from enum import Enum, auto
import copy
import pickle

import config

from env_var_names import *
from shell_ast.ast_util import *
from shasta.ast_node import ast_match, is_empty_cmd, string_of_arg
from shasta.json_to_ast import to_ast_node
from parse import from_ast_objects_to_shell
from speculative import util_spec

## There are two types of ast_to_ast transformations
class TransformationType(Enum):
    PASH = 'pash'
    SPECULATIVE = 'spec'

class ShellLoopContext:
    def __init__(self, test_block, next_block):
        self.test_block = test_block
        self.next_block = next_block

class ShellIfContext:
    def __init__(self, test_block, next_block, has_else=False):
        self.test_block = test_block
        self.next_block = next_block
        self.has_else = has_else

# class ShellWhileContext:
#     def __init__(self, test_block, next_block):
#         self.test_block = test_block
#         self.next_block = next_block

class ShellBB:
    commands: list[int] # not used atm
    def __init__(self, num: int):
        self.num = num
        self._is_emtpy = True
        # self.commands = []

    def add_command(self, command):
        self._is_emtpy = False

    def make_non_empty(self):
        self._is_emtpy = False

    def is_empty(self):
        return self._is_emtpy

class EdgeReason(Enum):
    IF_TAKEN = auto()
    ELSE_TAKEN = auto()
    LOOP_TAKEN = auto()
    LOOP_SKIP = auto()
    LOOP_BACK = auto()
    LOOP_BEGIN = auto()
    LOOP_END = auto()
    OTHER = auto()

    def __str__(self):
        return f'{self.name}'

class ShellProg:
    def __init__(self):
        self.ast_nodes = []
        self.init_bb = ShellBB(0)
        self.bbs = [self.init_bb]
        self.current_bb = 0
        self.edges = {}
        self.contexts = []

    def add_bb(self) -> int:
        next_bb = len(self.bbs)
        self.bbs.append(ShellBB(next_bb))
        return next_bb

    def add_edge(self, from_bb: int, to_bb: int, label: EdgeReason):
        assert from_bb >= 0 and from_bb < len(self.bbs)
        assert to_bb >= 0 and to_bb < len(self.bbs)
        if not from_bb in self.edges:
            self.edges[from_bb] = {}
        self.edges[from_bb][to_bb] = label

    def enter_for(self):
        test_bb = self.add_bb()
        self.add_edge(self.current_bb, test_bb, EdgeReason.LOOP_BEGIN)
        next_bb = self.add_bb()
        self.add_edge(test_bb, next_bb, EdgeReason.LOOP_SKIP)
        body_bb = self.add_bb()
        self.add_edge(test_bb, body_bb, EdgeReason.LOOP_TAKEN)
        self.contexts.append(ShellLoopContext(test_bb, next_bb))
        self.current_bb = body_bb

    def leave_for(self):
        loop_context = self.contexts.pop()
        assert isinstance(loop_context, ShellLoopContext)
        test_bb = loop_context.test_block
        next_bb = loop_context.next_block
        self.add_edge(self.current_bb, test_bb, EdgeReason.LOOP_BACK)
        self.current_bb = next_bb

    def enter_while(self):
        self.enter_for()

    def leave_while(self):
        self.leave_for()

    def enter_if(self):
        test_bb = self.current_bb
        next_bb = self.add_bb()
        body_bb = self.add_bb()
        self.add_edge(test_bb, body_bb, EdgeReason.IF_TAKEN)
        self.contexts.append(ShellIfContext(test_bb, next_bb))
        self.current_bb = body_bb

    def enter_else(self):
        if_context = self.contexts.pop()
        assert isinstance(if_context, ShellIfContext)
        test_bb = if_context.test_block
        next_bb = if_context.next_block
        self.add_edge(self.current_bb, next_bb, EdgeReason.OTHER)
        else_bb = self.add_bb()
        self.add_edge(test_bb, else_bb, EdgeReason.ELSE_TAKEN)
        self.contexts.append(ShellIfContext(test_bb, next_bb, has_else=True))
        self.current_bb = else_bb

    def leave_if(self):
        if_context = self.contexts.pop()
        assert isinstance(if_context, ShellIfContext)
        test_bb = if_context.test_block
        next_bb = if_context.next_block
        if not if_context.has_else:
            self.add_edge(test_bb, next_bb, EdgeReason.ELSE_TAKEN)
        self.add_edge(self.current_bb, next_bb, EdgeReason.OTHER)
        self.current_bb = next_bb

    def add_break(self):
        assert len(self.contexts) > 0
        for context in reversed(self.contexts):
            if isinstance(context, ShellLoopContext):
                break
        else:
            assert False
        assert isinstance(context, ShellLoopContext)
        break_bb = self.current_bb
        self.add_edge(break_bb, context.next_block, EdgeReason.LOOP_END)

    def add_command(self, command):
        self.bbs[self.current_bb].add_command(command)

## Use this object to pass options inside the preprocessing
## trasnformation.
class TransformationState:
    def __init__(self, mode: TransformationType):
        self.mode = mode
        self.node_counter = 0
        self.loop_counter = 0
        self.loop_contexts = []
        self.prog = ShellProg()

    def get_mode(self):
        return self.mode

    ## Node id related
    def get_next_id(self):
        new_id = self.node_counter
        self.node_counter += 1
        return new_id

    def get_current_id(self):
        return self.node_counter - 1

    def get_number_of_ids(self):
        return self.node_counter

    ## Loop id related
    def get_next_loop_id(self):
        new_id = self.loop_counter
        self.loop_counter += 1
        return new_id

    def get_current_loop_context(self):
        ## We want to copy that
        return self.loop_contexts[:]

    def get_current_loop_id(self):
        if len(self.loop_contexts) == 0:
            return None
        else:
            return self.loop_contexts[0]

    def current_bb(self):
        return self.prog.current_bb

    def enter_loop(self):
        new_loop_id = self.get_next_loop_id()
        self.loop_contexts.insert(0, new_loop_id)
        self.prog.enter_for()
        return new_loop_id

    def exit_loop(self):
        self.loop_contexts.pop(0)
        self.prog.leave_for()

    def enter_if(self):
        self.prog.enter_if()

    def enter_else(self):
        self.prog.enter_else()

    def exit_if(self):
        self.prog.leave_if()

    def visit_command(self, command):
        if len(command.arguments) > 0 and string_of_arg(command.arguments[0]) == 'break':
            self.prog.add_break()
        else:
            self.prog.add_command(command)

## TODO: Turn it into a Transformation State class, and make a subclass for
##       each of the two transformations. It is important for it to be state, because
##       it will need to be passed around while traversing the tree.
class SpeculativeTransformationState(TransformationState):
    def __init__(self, mode: TransformationType, po_file: str):
        super().__init__(mode)
        assert(self.mode is TransformationType.SPECULATIVE)
        self.partial_order_file = po_file
        self.partial_order_edges = []
        self.partial_order_node_bb = {}

    def get_partial_order_file(self):
        assert(self.mode is TransformationType.SPECULATIVE)
        return self.partial_order_file

    def add_edge(self, from_id: int, to_id: int):
        self.partial_order_edges.append((from_id, to_id))

    def get_all_edges(self):
        return self.partial_order_edges

    def add_node_bb(self, node_id: int, bb_id: int):
        self.partial_order_node_bb[node_id] = bb_id

    def get_all_node_bb(self):
        return self.partial_order_node_bb


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
        assert(isinstance(ast, AstNode))

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
        ## If the current candidate dataflow region is non-empty
        ## it means that the previous AST was in the background so
        ## the current one has to be included in the process no matter what
        elif (len(candidate_dataflow_region) > 0):
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
        elif(preprocessed_ast_object.should_replace_whole_ast()):
            replaced_ast = replace_df_region([preprocessed_ast_object.ast], trans_options,
                                             ast_text=original_text, disable_parallel_pipelines=last_object)
            preprocessed_asts.append(replaced_ast)

        ## In this case, it is possible that no replacement happened,
        ## meaning that we can simply return the original parsed text as it was.
        elif(preprocessed_ast_object.will_anything_be_replaced() or original_text is None):
            preprocessed_asts.append(preprocessed_ast_object.ast)
        elif trans_options.get_mode() is TransformationType.SPECULATIVE \
        and len(preprocessed_ast_object.ast.arguments) == 0 \
        and len(preprocessed_ast_object.ast.assignments) > 0:
            replaced_ast = replace_df_region([preprocessed_ast_object.ast], trans_options,
                                            ast_text=original_text, disable_parallel_pipelines=last_object)
            preprocessed_asts.append(replaced_ast)
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
    return ast_match(ast_object, preprocess_cases, trans_options, last_object)

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
    trans_options : TransformationState
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
    trans_options.visit_command(ast_node)
    log(f"[DEBUG_LOG] command: {ast_node}")
    return preprocessed_ast_object

# Background of (linno * t * redirection list)
## TODO: It might be possible to actually not close the inner node but rather apply the redirections on it
def preprocess_node_redir(ast_node, trans_options, last_object=False):
    preprocessed_node, something_replaced = preprocess_close_node(ast_node.node, trans_options,
                                                                  last_object=last_object)
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
    ## If we are in a loop, we push the loop identifier into the loop context
    loop_id = trans_options.enter_loop()
    preprocessed_body, something_replaced = preprocess_close_node(ast_node.body, trans_options, last_object=last_object)

    ## TODO: Then send this iteration identifier when talking to the spec scheduler
    ## TODO: After running checks put this behind a check to only run under speculation

    ## Create a new variable that tracks loop iterations
    var_name = loop_iter_var(loop_id)
    export_node = make_export_var_constant_string(var_name, '0')
    increment_node = make_increment_var(var_name)

    ## Also store the whole sequence of loop iters in a file
    all_loop_ids = trans_options.get_current_loop_context()

    ## export pash_loop_iters="$pash_loop_XXX_iter $pash_loop_YYY_iter ..."
    save_loop_iters_node = export_pash_loop_iters_for_current_context(all_loop_ids)

    ## Prepend the increment in the body
    ast_node.body = make_typed_semi_sequence(
        [to_ast_node(increment_node),
         to_ast_node(save_loop_iters_node),
         copy.deepcopy(preprocessed_body)])

    ## We pop the loop identifier from the loop context.
    ##
    ## KK 2023-04-27: Could this exit happen before the replacement leading to wrong
    ##     results? I think not because we use the _close_node preprocessing variant.
    ##     A similar issue might happen for while
    trans_options.exit_loop()

    ## reset the loop iters after we exit the loop
    out_of_loop_loop_ids = trans_options.get_current_loop_context()
    reset_loop_iters_node = export_pash_loop_iters_for_current_context(out_of_loop_loop_ids)

    ## Prepend the export in front of the loop
    # new_node = ast_node
    new_node = make_typed_semi_sequence(
        [to_ast_node(export_node),
         ast_node,
         to_ast_node(reset_loop_iters_node)])
    # print(new_node)

    preprocessed_ast_object = PreprocessedAST(new_node,
                                              replace_whole=False,
                                              non_maximal=False,
                                              something_replaced=something_replaced,
                                              last_ast=last_object)

    return preprocessed_ast_object

def preprocess_node_while(ast_node, trans_options, last_object=False):
    ## If we are in a loop, we push the loop identifier into the loop context
    trans_options.enter_loop()

    preprocessed_test, sth_replaced_test = preprocess_close_node(ast_node.test, trans_options, last_object=last_object)
    preprocessed_body, sth_replaced_body = preprocess_close_node(ast_node.body, trans_options, last_object=last_object)
    ast_node.test = preprocessed_test
    ast_node.body = preprocessed_body
    something_replaced = sth_replaced_test or sth_replaced_body
    preprocessed_ast_object = PreprocessedAST(ast_node,
                                              replace_whole=False,
                                              non_maximal=False,
                                              something_replaced=something_replaced,
                                              last_ast=last_object)

    ## We pop the loop identifier from the loop context.
    trans_options.exit_loop()
    return preprocessed_ast_object

## This is the same as the one for `For`
def preprocess_node_defun(ast_node, trans_options, last_object=False):
    ## TODO: For now we don't want to compile function bodies
    # preprocessed_body = preprocess_close_node(ast_node.body)
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
    ast_node.body = preprocessed_body
    preprocessed_ast_object = PreprocessedAST(ast_node,
                                              replace_whole=False,
                                              non_maximal=False,
                                              something_replaced=sth_replaced,
                                              last_ast=last_object)
    return preprocessed_ast_object


def preprocess_node_if(ast_node, trans_options, last_object=False):
    # preprocessed_left, should_replace_whole_ast, is_non_maximal = preprocess_node(ast_node.left, irFileGen, config)
    trans_options : TransformationState
    preprocessed_cond, sth_replaced_cond = preprocess_close_node(ast_node.cond, trans_options, last_object=last_object)
    trans_options.enter_if()
    preprocessed_then, sth_replaced_then = preprocess_close_node(ast_node.then_b, trans_options, last_object=last_object)
    if not is_empty_cmd(ast_node.else_b):
        trans_options.enter_else()
    preprocessed_else, sth_replaced_else = preprocess_close_node(ast_node.else_b, trans_options, last_object=last_object)
    trans_options.exit_if()
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
def replace_df_region(asts, trans_options, disable_parallel_pipelines=False, ast_text=None) -> AstNode:

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
        text_to_output = get_shell_from_ast(asts, ast_text=ast_text)
        ## Generate an ID
        df_region_id = trans_options.get_next_id()

        ## Get the current loop id and save it so that the runtime knows
        ## which loop it is in.
        loop_id = trans_options.get_current_loop_id()

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
        replaced_node = make_call_to_spec_runtime(df_region_id, loop_id)
    else:
        ## Unreachable
        assert(False)

    return to_ast_node(replaced_node)


def get_shell_from_ast(asts, ast_text=None) -> str:
    ## If we don't have the original ast text, we need to unparse the ast
    if (ast_text is None):
        text_to_output = from_ast_objects_to_shell(asts)
    else:
        text_to_output = ast_text
    return text_to_output


##
## Code that constructs the preprocessed ASTs
##


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

    ## Call the runtime
    arguments = [string_to_argument("source"),
                 string_to_argument(config.RUNTIME_EXECUTABLE)]
    runtime_node = make_command(arguments,
                                assignments=assignments)
    return runtime_node

## TODO: Make that an actual call to the spec runtime
def make_call_to_spec_runtime(command_id: int, loop_id) -> AstNode:
    assignments = [["pash_spec_command_id",
                        string_to_argument(str(command_id))]]
    if loop_id is None:
        loop_id_str = ""
    else:
        loop_id_str = str(loop_id)

    assignments.append(["pash_spec_loop_id",
                        string_to_argument(loop_id_str)])

    ## Call the runtime
    arguments = [string_to_argument("source"),
                 string_to_argument(config.RUNTIME_EXECUTABLE)]
    ## Pass all relevant argument to the planner
    runtime_node = make_command(arguments,
                                assignments=assignments)

    return runtime_node
