from definitions.ast_node_c import *
from definitions.no_match_exception import *

class AstNode:
    # create an AstNode object from an ast object as parsed by libdash
    def __init__(self, ast_object):
        try:
            self.construct = AstNodeConstructor(ast_object[0])
            self.parse_args(ast_object[1])
        except ValueError as no_matching_construct:
            raise NoMatchException('{} is not a construct we can handle'.format(ast_object[0]))

    def parse_args(self, args):
        if self.construct is AstNodeConstructor.PIPE:
            self.is_background = args[0]
            self.items = args[1]

        elif self.construct is AstNodeConstructor.COMMAND:
            self.line_number = args[0]
            self.assignments = args[1]
            self.arguments = args[2]
            self.redir_list = args[3]
            return

        elif self.construct is AstNodeConstructor.AND or AstNodeConstructor.OR or AstNodeConstructor.SEMI:
            self.left_operand = args[0]
            self.right_operand = args[1]

        elif self.construct is AstNodeConstructor.REDIR or AstNodeConstructor.SUBSHELL or AstNodeConstructor.BACKGROUND:
            self.line_number = args[0]
            # TODO maybe pick a better name?
            self.node = args[1]
            self.redir_list = args[2]

        elif self.construct is AstNodeConstructor.DEFUN:
            self.line_number = args[0]
            self.name = args[1]
            self.body = args[2]

    def check(self, **kwargs):
        # user-supplied custom checks
        for key, value in kwargs.items():
            try:
                assert(value())
            except Exception as exc:
                print("check for {} construct failed at key {}".format(self.construct, key))
                raise exc

