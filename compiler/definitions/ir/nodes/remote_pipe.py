from definitions.ir.dfg_node import *


class RemotePipe(DFGNode):
    def __init__(
        self,
        inputs,
        outputs,
        com_name,
        com_category,
        com_options=None,
        com_redirs=None,
        com_assignments=None,
    ):
        com_options = [] if com_options is None else com_options
        com_redirs = [] if com_redirs is None else com_redirs
        com_assignments = [] if com_assignments is None else com_assignments

        super().__init__(
            inputs,
            outputs,
            com_name,
            com_category,
            com_options=com_options,
            com_redirs=com_redirs,
            com_assignments=com_assignments,
        )


def make_remote_pipe(inputs, outputs, host_ip, port, is_remote_read, id):
    com_category = "pure"
    options = []
    opt_count = 0

    if is_remote_read:
        remote_pipe_bin = os.path.join(
            config.PASH_TOP, config.config["runtime"]["remote_read_binary"]
        )
    else:
        remote_pipe_bin = os.path.join(
            config.PASH_TOP, config.config["runtime"]["remote_write_binary"]
        )

    com_name = Arg.string_to_arg(remote_pipe_bin)

    options.append((opt_count, Arg.string_to_arg(f"--addr {host_ip}:{port}")))
    options.append((opt_count + 1, Arg.string_to_arg(f"--id {id}")))

    return RemotePipe(inputs, outputs, com_name, com_category, com_options=options)
