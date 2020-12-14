from definitions.ir.command import *
from definitions.ir.file_id import *

from ir_utils import string_to_argument, make_command

import config

import os

class Split(Command):
    def __init__(self, ast, command, options, in_stream, out_stream,
                 opt_indices, category, stdin=None, stdout=None, redirections=[]):
        super().__init__(ast, command, options, in_stream, out_stream,
                         opt_indices, category, stdin=stdin, stdout=stdout, redirections=redirections)
    
    def to_ast(self, _drain_streams):
        ## We don't care about the first index which is the batch size.
        opt_indices = self.all_opt_indices()[1:]
        option_asts = [self.options[opt_i].to_ast() for opt_i in opt_indices]
        ## Find the auto split binary in config
        auto_split_bin = os.path.join(config.PASH_TOP, config.config['runtime']['auto_split_binary'])
        ## Find the input fid
        ## TODO: Improve this. At the moment split actually takes its input file as a first argument
        ##       but in the IR we have it as stdin
        fid = Find(self.stdin)
        arguments = [string_to_argument(auto_split_bin)] + [fid.to_ast()] + option_asts
        node = make_command(arguments)
        return node

## TODO: Make a proper splitter subclass of Node
def make_split_file(in_fid, out_fids, batch_size):
    ## TODO: I probably have to give the file names as options to the command to.
    options = [string_to_argument(str(batch_size))] + out_fids
    opt_indices = [("option", 0)]
    stdout_indices = [("option", i+1) for i in range(len(out_fids))]
    command = Split(None, # TODO: Make a proper AST
                    string_to_argument("split_file"),
                    options,
                    ["stdin"],
                    stdout_indices,
                    opt_indices,
                    None, # TODO: Category?
                    in_fid)
    return command
