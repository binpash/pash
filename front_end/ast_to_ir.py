from ir import *
from union_find import *

import pickle

## Checks if the given ASTs are supported
def check_if_asts_supported(ast_objects):
    ## TODO: Implement
    return

##
## Compile AST -> Extended AST with IRs
##

## Compiles a given AST to an intermediate representation tree, which
## has some subtrees in it that are graphs representing a distributed
## computation.
##
## The above assumes that subtrees of the AST are disjoint
## computations that can be distributed separately (locally/modularly)
## without knowing about previous or later subtrees that can be
## distributed. Is that reasonable?
def compile_ast(ast_object, fileIdGen):
    compiled_ast = compile_node(ast_object, fileIdGen)
    return compiled_ast

def compile_node(ast_node, fileIdGen):
    # print("Compiling node: {}".format(ast_node))

    ## TODO: Can this become less verbose? (e.g. skip the lambda)
    cases = {
        "Pipe": (lambda c, b, i: compile_node_pipe(c, b, i, fileIdGen)),
        "Command": (lambda c, l, ass, args, r: compile_node_command(c, l, ass, args, r, fileIdGen)),
        "And" : (lambda c, l, r: compile_node_and_or_semi(c, l, r, fileIdGen)),
        "Or" : (lambda c, l, r: compile_node_and_or_semi(c, l, r, fileIdGen)),
        "Semi" : (lambda c, l, r: compile_node_and_or_semi(c, l, r, fileIdGen)),
        "Redir" : (lambda c, l, n, r: compile_node_redir(c, l, n, r, fileIdGen)),
        "Subshell" : (lambda c, l, n, r: compile_node_subshell(c, l, n, r, fileIdGen)),
        "Background" : (lambda c, l, n, r: compile_node_background(c, l, n, r, fileIdGen)),
        "Defun" : (lambda c, l, n, b: compile_node_defun(c, l, n, b, fileIdGen))
    }

    compiled_ast = ast_match(ast_node, cases)
    
    # print("Compiled node: {}".format(compiled_ast))
    return compiled_ast


def compile_node_pipe(construct, background, pipe_items, fileIdGen):
    ## Note: Background indicates when the pipe should be run in the background.
    ##
    ## TODO: Investigate whether we can optimize more by running
    ##       the background pipes in a distributed fashion.
    compiled_pipe_nodes = combine_pipe([compile_node(pipe_item, fileIdGen)
                                            for pipe_item in pipe_items])

    if (len(compiled_pipe_nodes) == 1):
        ## Note: When calling combine_pipe_nodes (which
        ##       optimistically distributes all the children of a
        ##       pipeline) the compiled_pipe_nodes should always
        ##       be one IR
        compiled_ir = compiled_pipe_nodes[0]
        ## Note: Save the old ast for the end-to-end prototype
        ast_node = {construct : [background, pipe_items]}
        compiled_ir.set_ast(ast_node)
        compiled_ast = compiled_ir
    else:
        ## Note: This should be unreachable with the current combine pipe
        assert(False)
        compiled_ast = {construct : [arguments[0]] + [compiled_pipe_nodes]}
    return compiled_ast

## This combines all the children of the Pipeline to an IR, even
## though they might not be IRs themselves. This means that an IR
## might contain mixed commands and ASTs. The ASTs can be
## (conservatively) considered as stateful commands by default).
def combine_pipe(ast_nodes):
    ## Initialize the IR with the first node in the Pipe
    if (isinstance(ast_nodes[0], IR)):
        combined_nodes = ast_nodes[0]
    else:
        combined_nodes = IR([ast_nodes[0]])

    ## Combine the rest of the nodes
    for ast_node in ast_nodes[1:]:
        if (isinstance(ast_node, IR)):
            combined_nodes.pipe_append(ast_node)
        else:
            ## FIXME: This one will not work. The IR of an AST node
            ##        doesn't have any stdin or stdout.
            combined_nodes.pipe_append(IR([ast_node]))
            
    return [combined_nodes]

def compile_node_command(construct, lineno, assignments, args, redir_list, fileIdGen):

    ast_node = {construct : [lineno, assignments, args, redir_list]}
    
    ## TODO: Do we need the line number?

    ## TODO: We probably also need to do something with the redirection list.
    
    compiled_assignments = compile_assignments(assignments, fileIdGen)
    
    ## If there are no arguments, the command is just an
    ## assignment
    if(len(args) == 0):
        ## Just compile the assignments. Specifically compile the
        ## assigned values, because they might have command
        ## substitutions etc..
        compiled_ast = {construct : [lineno] + [compiled_assignments] + [args, redir_list]}
    else:
        command_name = args[0]
        options = compile_command_arguments(args[1:], fileIdGen)
        
        stdin_fid = fileIdGen.next_file_id()
        stdout_fid = fileIdGen.next_file_id()
        ## Question: Should we return the command in an IR if one of
        ## its arguments is a command substitution? Meaning that we
        ## will have to wait for its command to execute first?
        compiled_ast = IR([Command(ast_node,
                                   command_name,
                                   options,
                                   stdin = stdin_fid,
                                   stdout = stdout_fid)],
                          stdin = stdin_fid,
                          stdout = stdout_fid)
        compiled_ast.set_ast(ast_node)
    
    return compiled_ast

def compile_node_and_or_semi(construct, left_node, right_node, fileIdGen):
    compiled_ast = {construct : [compile_node(left_node, fileIdGen),
                                 compile_node(right_node, fileIdGen)]}
    return compiled_ast

def compile_node_redir(construct, lineno, node, redir_list, fileIdGen):
    compiled_node = compile_node(node, fileIdGen)
    
    if (isinstance(compiled_node, IR)):
        ## TODO: I should use the redir list to redirect the files of
        ##       the IR accordingly
        compiled_ast = compiled_node
    else:
        compiled_ast = {construct : [lineno, compiled_node, redir_list]}

    return compiled_ast

def compile_node_subshell(construct, lineno, node, redir_list, fileIdGen):
    compiled_node = compile_node(node, fileIdGen)
    
    ## Question: It seems that subshell can be handled exactly
    ##           like a redir. Is that true?
    
    ## TODO: Make sure that propagating the IR up, doesn't create
    ##       any issue.
    
    if (isinstance(compiled_node, IR)):
        ## TODO: I should use the redir list to redirect the files of
        ##       the IR accordingly
        compiled_ast = compiled_node
    else:
        compiled_ast = {construct : [lineno, compiled_node, redir_list]}

    return compiled_ast

def compile_node_background(construct, lineno, node, redir_list, fileIdGen):
    compiled_node = compile_node(node, fileIdGen)
    
    ## TODO: I should use the redir list to redirect the files of
    ##       the IR accordingly
    if (isinstance(compiled_node, IR)):
        ## TODO: Redirect the stdout, stdin accordingly
        compiled_ast = compiled_node
    else:
        ## Note: It seems that background nodes can be added in
        ##       the distributed graph similarly to the children
        ##       of pipelines.
        ##
        ## Question: What happens with the stdin, stdout. Should
        ## they be closed?
        compiled_ast = IR([compiled_node])

    return compiled_ast

def compile_node_defun(construct, lineno, name, body, fileIdGen):
    ## It is not clear how we should handle functions.
    ##
    ## - Should we transform their body to IR?
    ## - Should we handle calls to the functions as commands?
    ##
    ## It seems that we should do both. But we have to think if
    ## this introduces any possible problem.
    
    ## TODO: Investigate whether it is fine to just compile the
    ##       body of functions.
    compiled_body = compile_node(body, fileIdGen)
    compiled_ast = {construct : [lineno, name, compiled_body]}
    return compiled_ast


def compile_arg_char(arg_char, fileIdGen):
    key, val = get_kv(arg_char)
    if (key == 'C'):
        return arg_char
    elif (key == 'B'):
        ## TODO: I probably have to redirect the input of the compiled
        ##       node (IR) to be closed, and the output to be
        ##       redirected to some file that we will use to write to
        ##       the command argument to complete the command
        ##       substitution.
        compiled_node = compile_node(val, fileIdGen)
        return {key : compiled_node}
    elif (key == 'Q'):
        compiled_val = compile_command_argument(val, fileIdGen)
        return {key : compiled_val}
    else:
        ## TODO: Complete this
        return arg_char
    
def compile_command_argument(argument, fileIdGen):
    compiled_argument = [compile_arg_char(char, fileIdGen) for char in argument]
    return compiled_argument
    
def compile_command_arguments(arguments, fileIdGen):
    compiled_arguments = [compile_command_argument(arg, fileIdGen) for arg in arguments]
    return compiled_arguments

## Compiles the value assigned to a variable using the command argument rules.
## TODO: Is that the correct way to handle them?
def compile_assignments(assignments, fileIdGen):
    compiled_assignments = [[assignment[0], compile_command_argument(assignment[1], fileIdGen)]
                            for assignment in assignments]
    return compiled_assignments




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
##
## TODO: Optimization: This pass could be merged with the previous
## translation one. I am not doing this now so that it remains
## cleaner. If we find out that we might benefit from this
## optimization, we can do it.
##
## TODO: Visual improvement: Can we abstract away the traverse of this
## ast in a mpa function, so that we don't rewrite all the cases?
def replace_irs(ast, irFileGen):

    if (isinstance(ast, IR)):
        replaced_ast = replace_ir(ast, irFileGen)
    else:
    
        cases = {
            ## Note: We should never encounter a Pipe construct, since all
            ## of them must have been become IRs
            ##
            ## "Pipe": (lambda c, b, i: {c : [b, i]}),
            "Command": (lambda c, l, ass, args, r: replace_irs_command(c, l, ass, args, r, irFileGen)),
            
            "And" : (lambda c, l, r: replace_irs_and_or_semi(c, l, r, irFileGen)),
            "Or" : (lambda c, l, r: replace_irs_and_or_semi(c, l, r, irFileGen)),
            "Semi" : (lambda c, l, r: replace_irs_and_or_semi(c, l, r, irFileGen)),
            
            ## TODO: Complete these
            # "Redir" : (lambda c, l, n, r: compile_node_redir(c, l, n, r, fileIdGen)),
            # "Subshell" : (lambda c, l, n, r: compile_node_subshell(c, l, n, r, fileIdGen)),
            # "Background" : (lambda c, l, n, r: compile_node_background(c, l, n, r, fileIdGen)),
            # "Defun" : (lambda c, l, n, b: compile_node_defun(c, l, n, b, fileIdGen))
        }
        
        replaced_ast = ast_match(ast, cases)
    
    return replaced_ast

## This function serializes an IR in a file, and in its place, it adds
## a command that calls our distribution planner with the name of the
## saved file.
def replace_ir(ast_node, irFileGen):
    ir_file_id = irFileGen.next_file_id()

    ## TODO: Remember to not have the prefix hardocded :( In order to
    ## do this properly, I have to set this directory up when
    ## execution starts.

    ## TODO: I probably have to generate random file names for the irs, so
    ## that multiple executions of dish don't interfere.
    ir_filename = ir_file_id.toFileName("/tmp/dish_temp_ir")
    
    ## Serialize the node in a file
    with open(ir_filename, "wb") as ir_file:
        pickle.dump(ast_node, ir_file)
    
    ## Replace it with a command that calls the distribution
    ## planner with the name of the file.
    replaced_node = make_command(ir_filename)
    
    return replaced_node

## This function makes a command that calls the distribution planner
## together with the name of the file containing an IR. Then the
## distribution planner should read from this file and continue
## execution.
##
## TODO: Make sure stuff are not hardcoded in here
##
## (MAYBE) TODO: The way I did it, is by calling the parser once, and seeing
## what it returns. Maybe it would make sense to call the parser on
## the fly to have a cleaner implementation here?
def make_command(ir_filename):

    ## TODO: Do we need to do anything with the lineno? If so, make
    ## sure that I keep it in the IR, so that I can find it.
    arguments = [string_to_argument("python3"),
                 string_to_argument("distr_plan.py"),
                 string_to_argument(ir_filename)]
    lineno = 0
    node = { 'Command' : [lineno, [], arguments, []] }
    return node
                 

def replace_irs_and_or_semi(c, left, right, irFileGen):
    r_left = replace_irs(left, irFileGen)
    r_right = replace_irs(right, irFileGen)
    return {c : [r_left, r_right]}

def replace_irs_command(c, lineno, assignments, args, redir_list, irFileGen):
    ## TODO: I probably have to also handle the redir_list (applies to
    ## compile_node_command too)

    replaced_assignments = replace_irs_assignments(assignments, irFileGen)
    if(len(args) == 0):
        ## Note: The flow here is similar to compile
        replaced_ast = {c : [lineno, replaced_assignments, args, redir_list]}
    else:
        command_name = args[0]
        replaced_options = replace_irs_command_arguments(args[1:], irFileGen)
        replaced_args = [command_name] + replaced_options
        
        replaced_ast = {c : [lineno, replaced_assignments, replaced_args, redir_list]}

    return replaced_ast

def replace_irs_arg_char(arg_char, irFileGen):
    key, val = get_kv(arg_char)
    if (key == 'C'):
        return arg_char
    elif (key == 'B'):
        replaced_node = replace_irs(val, irFileGen)
        return {key : replaced_node}
    elif (key == 'Q'):
        replaced_val = replace_irs_command_argument(val, irFileGen)
        return {key : replaced_val}
    else:
        ## TODO: Complete this (as we have to do with the compile)
        return arg_char
    
def replace_irs_command_argument(argument, irFileGen):
    replaced_argument = [replace_irs_arg_char(char, irFileGen) for char in argument]
    return replaced_argument
    
def replace_irs_command_arguments(arguments, irFileGen):
    replaced_arguments = [replace_irs_command_argument(arg, irFileGen) for arg in arguments]
    return replaced_arguments

## This is similar to compile_assignments
##
## TODO: Is that the correct way to handle them? Question also applied
## to compile_assignments.
def replace_irs_assignments(assignments, irFileGen):
    replaced_assignments = [[assignment[0], replace_irs_command_argument(assignment[1], irFileGen)]
                            for assignment in assignments]
    return replaced_assignments


    
##
## Pattern matching for the AST
##


## For now these checks are too simple. 
##
## Maybe we can move them to the check_if_ast_is_supported?
def check_pipe(construct, arguments):
    assert(len(arguments) == 2)
    ## The pipe should have at least 2 children
    assert(len(arguments[1]) >= 2)

def check_command(construct, arguments):
    assert(len(arguments) == 4)

def check_and(construct, arguments):
    assert(len(arguments) == 2)

def check_or(construct, arguments):
    assert(len(arguments) == 2)

def check_semi(construct, arguments):
    assert(len(arguments) == 2)

def check_redir(construct, arguments):
    assert(len(arguments) == 3)

def check_subshell(construct, arguments):
    assert(len(arguments) == 3)

def check_background(construct, arguments):
    assert(len(arguments) == 3)

def check_defun(construct, arguments):
    assert(len(arguments) == 3)

def get_case(constructor, cases):
    ##TODO: Throw a specific error, that indicates non-exhaustive
    ##pattern match
    return cases[constructor]

## This function does the pattern match. It throws an error if the
## pattern match is non-exhaustive, and it checks that the ast is
## indeed well-formed.
##
## TODO: Throw an error that says that pattern match failed in
## arguments, if arguments where malformed. To do this I should have
## everything in a try-catch assertion errors. Actually I can probably
## just inline all the check commands in the ast match function.
def ast_match(ast_node, cases):
    # print("Compiling node: {}".format(ast_node))

    construct, arguments = get_kv(ast_node)
    case = get_case(construct, cases)
    
    if (construct == 'Pipe'):
        check_pipe(construct, arguments)

        background = arguments[0]
        pipe_items = arguments[1]
        result = case(construct, background, pipe_items)

    elif (construct == 'Command'):
        check_command(construct, arguments)
        
        lineno = arguments[0]
        assignments = arguments[1]
        args = arguments[2]
        redir_list = arguments[3]
        result = case(construct, lineno, assignments, args, redir_list)
        
    elif (construct == 'And'):
        check_and(construct, arguments)
        
        left_node = arguments[0]
        right_node = arguments[1]
        result = case(construct, left_node, right_node)

    elif (construct == 'Or'):
        check_or(construct, arguments)
        
        left_node = arguments[0]
        right_node = arguments[1]
        result = case(construct, left_node, right_node)
        
    elif (construct == 'Semi'):
        check_semi(construct, arguments)
        
        left_node = arguments[0]
        right_node = arguments[1]
        result = case(construct, left_node, right_node)

    elif (construct == 'Redir'):
        check_redir(construct, arguments)

        line_no = arguments[0]
        node = arguments[1]
        redir_list = arguments[2]
        result = case(construct, line_no, node, redir_list)

    elif (construct == 'Subshell'):
        check_subshell(construct, arguments)

        line_no = arguments[0]
        node = arguments[1]
        redir_list = arguments[2]
        result = case(construct, line_no, node, redir_list)
            
    elif (construct == 'Background'):
        check_background(construct, arguments)

        line_no = arguments[0]
        node = arguments[1]
        redir_list = arguments[2]
        result = case(construct, line_no, node, redir_list)

    elif (construct == 'Defun'):
        check_defun(construct, arguments)
        
        line_no = arguments[0]
        name = arguments[1]
        body = arguments[2]
        result = case(construct, line_no, name, body)
        
    else:
        raise TypeError("Unimplemented construct: {}".format(construct))

    # print("Compiled node: {}".format(compiled_ast))
    return result
