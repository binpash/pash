from definitions.ir.node import *
from definitions.ir.arg import *
from definitions.ir.file_id import *
from definitions.ir.redirection import *

from ir_utils import *
from util import *
import config

import os

## Commands are specific Nodes that can be parallelized if they are
## classified as stateless, etc...
class Command(Node):
    def __init__(self, ast, command, options, in_stream, out_stream,
                 opt_indices, category, stdin=None, stdout=None, redirections=[]):
        super().__init__(ast, in_stream, out_stream, category, stdin, stdout)
        self.command = Arg(command)
        self.options = [opt if isinstance(opt, FileId) else Arg(opt) for opt in options]
        self.opt_indices = opt_indices
        self.redirections = [Redirection(redirection) for redirection in redirections]
        self.apply_redirections()

    def __repr__(self):
        prefix = "Command"
        if (self.category == "stateless"):
            prefix = "Stateless"
        elif (self.category == "pure"):
            prefix = "Pure"
        # output = "{}: \"{}\" in:{} out:{} opts:{}".format(
        #     prefix, self.command, self.stdin, self.stdout, self.options)
        output = "{}: \"{}\" in:{} out:{}".format(
            prefix, self.command, self.get_input_file_ids(),
            self.get_output_file_ids())
        return output

    def serialize(self):
        all_opt_indices = self.all_opt_indices()
        options_string = " ".join([self.options[opt_i].opt_serialize() for opt_i in all_opt_indices])
        output = "{} {}".format(self.command, options_string)
        return output


    ## TODO: Improve this functio to be separately implemented for different special nodes,
    ##       such as cat, eager, split, etc...
    def to_ast(self, drain_streams):    
        ## TODO: We might not want to implement this at all actually
        if (drain_streams):
            raise NotImplementedError()
        else:
            redirs = []
            all_opt_indices = self.all_opt_indices()
            option_asts = [self.options[opt_i].to_ast() for opt_i in all_opt_indices]
            # log("Options:", option_asts)
            arguments = [self.command.to_ast()] + option_asts 

            ## If the command has stdin, redirect the pipe to stdin                
            if ("stdin" in self.in_stream):
                fid = Find(self.stdin)
                redirs.append(redir_file_to_stdin(fid.to_ast()))

            ## If the command has stdout, redirect stdout to a pipe                
            if ("stdout" in self.out_stream):
                fid = Find(self.stdout)
                redirs.append(redir_stdout_to_file(fid.to_ast()))

            node = make_command(arguments,redirections=redirs)
            
        return node

    
    ## Gets all option indices
    def all_opt_indices(self):
        all_opt_indices = [o_i[1] for o_i in (self.opt_indices + self.in_stream + self.out_stream)
                           if isinstance(o_i, tuple)]
        all_opt_indices.sort()
        return all_opt_indices


    def get_non_file_options(self):
        return [self.options[i] for _, i in self.opt_indices]


    ## TODO: Maybe delete
    def get_file_id(self, chunk):
        if (isinstance(chunk, tuple)
                and len(chunk) == 2
                and chunk[0] == 'option'):
            return self.options[chunk[1]]

        return super().get_file_id(chunk)

    def set_file_id(self, chunk, value):
        if (isinstance(chunk, tuple)
                and len(chunk) == 2
                and chunk[0] == 'option'):
            self.options[chunk[1]] = value
            return

        super().set_file_id(chunk, value)
