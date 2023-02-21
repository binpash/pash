from enum import Enum

class AstNodeConstructor(Enum):
    PIPE = 'Pipe'
    COMMAND = 'Command'
    AND = 'And'
    OR = 'Or'
    NOT = 'Not'
    SEMI = 'Semi'
    REDIR = 'Redir'
    SUBSHELL = 'Subshell'
    BACKGROUND = 'Background'
    DEFUN = 'Defun'
    FOR = 'For'
    WHILE = 'While'
    IF = 'If'
    CASE = 'Case'