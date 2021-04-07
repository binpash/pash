from definitions.ir.arg import *
from ir_utils import *

class Redirection():
    def __init__(self, redirection):
        ## Handle initialization from an existing Redirection object
        if(isinstance(redirection, Redirection)):
            redirection = redirection.to_ast()
        assert(len(redirection) == 2)
        self.redir_type = redirection[0]
        assert(len(redirection[1]) == 3)
        self.redir_subtype = redirection[1][0]
        self.stream_id = redirection[1][1]
        self.file_arg = Arg(redirection[1][2])

        # log(redirection)
        ## TODO: Support all redirections
        assert(self.redir_type == 'File')
        assert(self.redir_subtype in ['To', 'From'])

    def __repr__(self):
        return '({}, {}, {}, {})'.format(self.redir_type,
                                         self.redir_subtype,
                                         self.stream_id,
                                         self.file_arg)

    def to_ast(self):
        redir = make_kv(self.redir_type,
                        [self.redir_subtype,
                         self.stream_id,
                         self.file_arg.to_ast()])
        return redir

    def is_to_file(self):
        return (self.redir_type == 'File'
                and self.redir_subtype == 'To')

    def is_for_stdout(self):
        return (self.stream_id == 1)

    def is_from_file(self):
        return (self.redir_type == 'File'
                and self.redir_subtype == 'From')

    def is_for_stdin(self):
        return (self.stream_id == 0)
