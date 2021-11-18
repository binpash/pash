from definitions.ast_node import *

## This class is used by the preprocessor in ast_to_ir
class PreprocessedAST:
    def __init__(self, ast, replace_whole, non_maximal, something_replaced=True, last_ast=False):
        assert(isinstance(ast, AstNode))
        self.ast = ast
        self.replace_whole = replace_whole
        self.non_maximal = non_maximal
        self.something_replaced = something_replaced
        self.last_ast = last_ast

    def should_replace_whole_ast(self):
        return self.replace_whole

    def is_non_maximal(self):
        return self.non_maximal
    
    def will_anything_be_replaced(self):
        return self.something_replaced

    def is_last_ast(self):
        return self.last_ast

## This class represents text that was not modified at all by preprocessing, and therefore does not
## need to be unparsed.
class UnparsedScript:
    def __init__(self, text):
        self.text = text