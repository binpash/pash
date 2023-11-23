from abc import ABC, abstractmethod
import os
from datetime import datetime

from util import *
from shell_ast.ast_util import *
from parse import from_ast_objects_to_shell
import config


class ExportOption(ABC):
    @abstractmethod
    def export_ast(self, asts):
        pass

    @abstractmethod
    def export_setup(self, ir, args):
        pass

    @abstractmethod
    def export_cleanup(self, ir, args):
        pass


class PaShExport(ExportOption):
    RM = "rm_pash_fifos"
    MKFIFO = "mkfifo_pash_fifos"

    def export_ast(self, asts):
        return from_ast_objects_to_shell(asts)

    def export_setup(self, ir, args):
        ephemeral_fids = [fid for fid in ir.all_fids() if fid.is_ephemeral()]

        asts = []
        ## Create an `rm -f` for each ephemeral fid
        rm_asts = [self._make_rm_f_ast(fid.to_ast()) for fid in ephemeral_fids]
        defun_rm_pash_fifos = make_defun(self.RM, make_semi_sequence(rm_asts))

        ## Create a `mkfifo` for each ephemeral fid
        mkfifo_asts = [self._make_mkfifo_ast(fid.to_ast()) for fid in ephemeral_fids]
        defun_mkfifos = make_defun(self.MKFIFO, make_semi_sequence(mkfifo_asts))

        call_rm_pash_fifos = make_command([string_to_argument(self.RM)])
        call_mkfifos = make_command([string_to_argument(self.MKFIFO)])

        asts = [defun_rm_pash_fifos, defun_mkfifos, call_rm_pash_fifos, call_mkfifos]
        return [to_ast_node(ast) for ast in asts]

    def export_cleanup(self, ir, args):
        clean_up_graph = args.termination == "clean_up_graph"
        log_file = args.log_file

        asts = []
        if clean_up_graph:
            ## TODO: Wait for all output nodes not just one
            pids = [[standard_var_ast("!")]]
            clean_up_path_script = os.path.join(
                config.PASH_TOP, config.config["runtime"]["clean_up_graph_binary"]
            )
            com_args = [
                string_to_argument("source"),
                string_to_argument(clean_up_path_script),
            ] + pids
            if log_file == "":
                com = make_command(com_args)
            else:
                redirection = redir_append_stderr_to_string_file(log_file)
                com = make_command(com_args, redirections=[redirection])
            asts = [com]
        else:
            ## Otherwise we just wait for all processes to die.
            wait_com = make_command([string_to_argument("wait")])
            exit_status = make_command([string_to_argument("internal_exec_status=$?")])
            asts = [wait_com, exit_status]

        ## Create an `rm -f` for each ephemeral fid
        call_rm_pash_funs = make_command([string_to_argument(self.RM)])
        asts.append(call_rm_pash_funs)

        ## Make the following command:
        #    (exit $internal_exec_status)
        exit_ec_ast = self._make_exit_ec_ast()
        asts.append(exit_ec_ast)

        return [to_ast_node(ast) for ast in asts]

    def _make_exit_ec_ast(
        self,
    ):
        command = make_command(
            [string_to_argument("exit"), [make_quoted_variable("internal_exec_status")]]
        )
        return make_subshell(command)

    def _make_rm_f_ast(self, arg):
        return make_command([string_to_argument("rm"), string_to_argument("-f"), arg])

    def _make_mkfifo_ast(self, arg):
        return make_command([string_to_argument("mkfifo"), arg])


class AirflowExport(ExportOption):
    def export_ast(self, asts):
        shell_list = []
        for ast in asts:
            if isinstance(ast, UnparsedScript):
                shell_list.append(
                    f"command_task = BashOperator(task_id='command_task', bash_command='{ast.text}'"
                )
            elif isinstance(ast, IfNode):
                shell_list.append(
                    f"cond_task = BashOperator(task_id='cond_task', bash_command='{ast.cond}), xcom_push=True\n"
                    """
                        @task.branch(task_id='branch')
                        def branch_func(ti=None):
                        xcom_value = bool(ti.xcom_pull(task_ids='cond_task'))
                        if xcom_value:
                           return 'then_task'
                        else:
                           return 'else_task'
                        """
                    f"then_task = BashOperator(task_id='then_task', bash_command='{ast.then_b}'"
                    f"else_task = BashOperator(task_id='then_task', bash_command='{ast.else_b}'"
                )
            else:
                shell_list.append(ast.pretty())

        return "\n".join(shell_list) + "\n"

    def export_setup(self, ir, args):
        return [make_command([string_to_argument("echo setup")])]

    def export_cleanup(self, ir, args):
        return [make_command([string_to_argument("echo cleanup")])]


def to_shell(ir, args):
    backend_start_time = datetime.now()

    drain_streams = args.termination == "drain_stream"
    export_opt = AirflowExport()

    # First call an IR to AST compilation pass
    # TODO: If we have just a single node maybe we should just instantiate that without anything else.
    # This is implemented below, but there is no need to turn it on now.
    #
    # If we only have a single node, then we don't need to make prologues, epilogues, etc
    # if len(ir.nodes) == 1:
    #     nodes = list(ir.nodes.values())
    #     assert(len(nodes) == 1)
    #     node = nodes[0]
    #     body = [node.to_ast(ir.edges, drain_streams)]
    #     return body

    # NOTE: We first need to make the main body because it might create additional ephemeral fids.
    body = ir.to_ast(drain_streams)
    setup = export_opt.export_setup(ir, args)
    cleanup = export_opt.export_cleanup(ir, args)
    output_asts = setup + body + cleanup

    # Then just call the parser.
    output_script = export_opt.export_ast(output_asts)
    print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    print(output_script)
    print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++")

    backend_end_time = datetime.now()
    print_time_delta("Backend", backend_start_time, backend_end_time)

    return output_script
