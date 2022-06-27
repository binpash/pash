from definitions.ir.dfg_node import *
from ir_utils import *

class Eager(DFGNode):
    def __init__(self, inputs, outputs, com_name, com_category, com_options = [], 
                 com_redirs = [], com_assignments=[],
                 intermediate = None):
        # BEGIN ANNO : hack for intermediate at the end
        self.intermediate = intermediate
        # END ANNO
        super().__init__(inputs, outputs, com_name, com_category,
                         com_options=com_options, 
                         com_redirs=com_redirs, 
                         com_assignments=com_assignments)

    # BEGIN ANNO : copied from DFG node for hack for intermediate at the end
    def to_ast(self, edges, drain_streams):
        log(f'do we get here?')
        ## TODO: We might not want to implement this at all actually
        if (drain_streams):
            raise NotImplementedError()
        else:
            ## Handle special node to ast here
            # node = self.special_to_ast(edges)
            # if node is not None:
            #     return node

            redirs = self._to_ast_aux_get_redirs()
            assignments = self.com_assignments
            ## Start filling in the arguments
            opt_arguments = []
            # BEGIN ANNO
            # get_command_invocation_prefix_from_dfg_node
            log(f'com_name: {self.com_name}')
            log(f'edges: {edges}')
            log(f'inputs: {self.inputs}')
            log(f'outputs: {self.outputs}')
            log(f'com_redirs: {self.com_redirs}')
            log(f'pos config: {self.positional_config_list}')
            log(f'pos input: {self.positional_input_list}')
            log(f'pos output: {self.positional_output_list}')
            log(f'com_options: {self.com_options}')
            log(f'flag_option_list: {self.flag_option_list}')

            # if self.implicit_use_of_stdin: # need to recompute
            # cat a list of inputs into it; redirect a single one
            # else:

            # OLD
            # for i, opt in self.com_options:
            #     ## Pad the argument list with None
            #     opt_arguments = pad(opt_arguments, i)
            #     opt_arguments[i] = opt.to_ast()
            # log(f'opt_arguments: {format_args([val for val in opt_arguments if val is not None])}')
            # NEW
            opt_arguments_new = [get_ast_for_flagoption(flagoption) for flagoption in self.flag_option_list]
            opt_arguments_new += [get_ast_for_argstringtype(arg) for arg in self.positional_config_list]
            log(f'opt_arguments_new: {format_args(opt_arguments_new)}')
            # END ANNO

            com_name_ast = self.com_name.to_ast()
            option_asts = [opt.to_ast() for _, opt in self.com_options]

            ##
            ## 1. Find the input and output fids
            ## 2. Construct the rest of the arguments and input/output redirections according to
            ##    the command IO
            input_fids = [edges[in_id][0] for in_id in self.get_input_list()]
            output_fids = [edges[out_id][0] for out_id in self.outputs]
            rest_argument_fids, new_redirs = create_command_arguments_redirs(com_name_ast,
                                                                             option_asts,
                                                                             input_fids,
                                                                             output_fids)

            ## Transform the rest of the argument fids to arguments
            ## Since some of the rest_arguments can be None (they only contain inputs and outputs)
            ## we need to make sure that we don't turn None objects to asts.
            ##
            ## The None fields need to be filtered out because they are taken care of by the interleave function.
            ##
            ## TODO: Is this actually OK?
            rest_arguments = [fid.to_ast()
                              for fid in rest_argument_fids
                              if not fid is None]
            log(f'rest_arguments: {format_args(rest_arguments)}')

            ## Interleave the arguments since options args might contain gaps.
            # BEGIN ANNO
            rest_arguments_backup = rest_arguments.copy()
            # OLD
            # arguments = interleave_args(opt_arguments, rest_arguments)
            # log(f'arguments fin: {format_args(arguments)}')
            # NEW
            arguments_new = opt_arguments_new + rest_arguments_backup + [self.intermediate.to_ast()]
            log(f'arguments_new: {format_args(arguments_new)}')
            # END ANNO

            all_arguments = [com_name_ast] + arguments_new
            all_redirs = redirs + new_redirs

            node = make_command(all_arguments, redirections=all_redirs, assignments=assignments)
        return node


def make_eager_node(input_id, output_id, intermediate_file_id, eager_exec_path):
    com_name = Arg(string_to_argument(eager_exec_path))
    com_category = "pure"
    ## TODO: In theory the intermediate file id is also an output...
    # BEGIN ANNO
    # OLD
    intermediate_identifier = Arg(intermediate_file_id.to_ast())
    com_options = [(2, intermediate_identifier)]
    return Eager([input_id],
                 [output_id],
                 com_name, 
                 com_category,
                 com_options=com_options,
                 intermediate=intermediate_identifier)
