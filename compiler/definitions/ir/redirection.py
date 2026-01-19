from definitions.ir.arg import *
from shell_ast.ast_util import *
from custom_error import UnparallelizableError

class Redirection:
    def __init__(self, redirection: RedirectionNode):
        ## TODO: Support all redirections
        if not isinstance(redirection, FileRedirNode) or redirection.redir_type not in [
            "To",
            "From",
        ]:
            raise NotImplementedError(redirection)

        self.redir_type = redirection.NodeName
        self.redir_subtype = redirection.redir_type
        self.stream_id = redirection.fd
        self.file_arg = Arg(redirection.arg)

        # log(redirection)

    def __repr__(self):
        return "({}, {}, {}, {})".format(
            self.redir_type, self.redir_subtype, self.stream_id, self.file_arg
        )

    def to_ast(self):
        redir = make_kv(
            self.redir_type,
            [self.redir_subtype, self.stream_id, self.file_arg.to_ast()],
        )
        return redir

    def is_to_file(self):
        return self.redir_type == "File" and self.redir_subtype == "To"

    def is_for_stdout(self):
        return self.stream_id == 1

    def is_from_file(self):
        return self.redir_type == "File" and self.redir_subtype == "From"

    def is_for_stdin(self):
        return self.stream_id == 0
