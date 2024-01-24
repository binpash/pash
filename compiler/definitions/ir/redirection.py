from definitions.ir.arg import *
from shell_ast.ast_util import *


class Redirection:
    def __init__(self, redirection: RedirectionNode):
        if isinstance(redirection, FileRedirNode):
            self.redir_type = FileRedirNode.NodeName
        elif isinstance(redirection, DupRedirNode):
            self.redir_type = DupRedirNode.NodeName
        elif isinstance(redirection, HeredocRedirNode):
            self.redir_type = HeredocRedirNode.NodeName

        self.redir_subtype = redirection.redir_type
        self.stream_id = redirection.fd
        self.file_arg = Arg(redirection.arg)

        # log(redirection)
        ## TODO: Support all redirections
        assert self.redir_type == "File"
        assert self.redir_subtype in ["To", "From"]

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
