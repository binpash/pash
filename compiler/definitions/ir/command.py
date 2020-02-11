from definitions.ir.node import *
from definitions.ir.arg import *
from definitions.ir.file_id import *

from ir_utils import *

## Commands are specific Nodes that can be parallelized if they are
## classified as stateless, etc...
class Command(Node):
    def __init__(self, ast, command, options, in_stream, out_stream,
                 opt_indices, category, stdin=None, stdout=None):
        super().__init__(ast, in_stream, out_stream, category, stdin, stdout)
        self.command = Arg(command)
        self.options = [opt if isinstance(opt, FileId) else Arg(opt) for opt in options]
        self.opt_indices = opt_indices

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
        all_opt_indices = [o_i[1] for o_i in (self.opt_indices + self.in_stream + self.out_stream)
                           if isinstance(o_i, tuple)]
        all_opt_indices.sort()
        options_string = " ".join([self.options[opt_i].opt_serialize() for opt_i in all_opt_indices])
        output = "{} {}".format(self.command, options_string)
        return output

    def get_non_file_options(self):
        return [self.options[i] for _, i in self.opt_indices]

    ## Get the file names of the outputs of the map commands. This
    ## differs if the command is stateless, pure that can be
    ## written as a map and a reduce, and a pure that can be
    ## written as a generalized map and reduce.
    def get_map_output_files(self, input_file_ids, fileIdGen):
        assert(self.category == "stateless" or self.is_pure_parallelizable())
        if(self.category == "stateless"):
            return [fileIdGen.next_file_id() for in_fid in input_file_ids]
        elif(self.is_pure_parallelizable()):
            return self.pure_get_map_output_files(input_file_ids, fileIdGen)
        else:
            print("Unreachable code reached :(")
            assert(False)
            ## This should be unreachable

    def pure_get_map_output_files(self, input_file_ids, fileIdGen):
        assert(self.is_pure_parallelizable())
        if(str(self.command) == "sort"):
            new_output_file_ids = [[fileIdGen.next_file_id()] for in_fid in input_file_ids]
        elif(str(self.command) == "bigrams_aux"):
            new_output_file_ids = [[fileIdGen.next_file_id() for i in range(BigramGMap.num_outputs)]
                                   for in_fid in input_file_ids]
        elif(str(self.command) == "alt_bigrams_aux"):
            new_output_file_ids = [[fileIdGen.next_file_id()] for in_fid in input_file_ids]
        else:
            print("Unreachable code reached :(")
            assert(False)
            ## This should be unreachable
        return new_output_file_ids

    def duplicate(self, old_input_file_id, new_input_file_ids, new_output_file_ids, fileIdGen):
        assert(self.category == "stateless" or self.is_pure_parallelizable())
        if(self.category == "stateless"):
            return self.stateless_duplicate(old_input_file_id, new_input_file_ids, new_output_file_ids)
        elif(self.is_pure_parallelizable()):
            return self.pure_duplicate(old_input_file_id, new_input_file_ids, new_output_file_ids, fileIdGen)
        else:
            print("Unreachable code reached :(")
            assert(False)
            ## This should be unreachable

    def stateless_duplicate(self, old_input_file_id, input_file_ids, output_file_ids):
        assert(self.category == "stateless")

        out_edge_file_ids = self.get_output_file_ids()
        assert(len(out_edge_file_ids) == 1)
        out_edge_file_id = out_edge_file_ids[0]

        ## Make a new cat command to add after the current command.
        new_cat = make_cat_node(output_file_ids, out_edge_file_id)

        in_out_file_ids = zip(input_file_ids, output_file_ids)

        new_commands = [self.find_in_out_and_make_duplicate_command(old_input_file_id, in_fid,
            out_edge_file_id, out_fid) for in_fid, out_fid in in_out_file_ids]

        return new_commands + [new_cat]

    def is_pure_parallelizable(self):
        return (self.category == "pure" and str(self.command) in list(map(get_command_from_definition, parallelizable_pure_commands)))

    def pure_duplicate(self, old_input_file_id, input_file_ids, output_file_ids, fileIdGen):
        assert(self.is_pure_parallelizable())

        in_out_file_ids = zip(input_file_ids, output_file_ids)

        simple_map_pure_commands = ["sort",
                                    "alt_bigrams_aux"]
        ## This is the category of all commands that don't need a
        ## special generalized map
        if(str(self.command) in simple_map_pure_commands):

            ## make_duplicate_command duplicates a node based on its
            ## output file ids, so we first need to assign them.
            new_output_file_ids = [fids[0] for fids in output_file_ids]

            new_commands = [self.find_in_out_and_make_duplicate_command(old_input_file_id, in_fid,
                self.stdout, out_fids[0]) for in_fid, out_fids in in_out_file_ids]
        elif(str(self.command) == "bigrams_aux"):
            new_commands = [BigramGMap([in_fid] + out_fids)
                            for in_fid, out_fids in in_out_file_ids]
        else:
            print("Unreachable code reached :(")
            assert(False)
            ## This should be unreachable
        return new_commands


    ## TODO: This has to change, to not search for the inputs in the
    ## in_stream
    def find_in_out_and_make_duplicate_command(self, old_input_file_id, in_fid,
            old_output_file_id, out_fid):

        ## First find what does the new file identifier refer to
        ## (stdin, or some argument)
        new_in_stream_index = self.find_file_id_in_in_stream(old_input_file_id)
        new_out_stream_index = self.find_file_id_in_out_stream(old_output_file_id)
        new_input_location = self.in_stream[new_in_stream_index]
        new_output_location = self.out_stream[new_out_stream_index]
        return self.make_duplicate_command(in_fid, out_fid, new_input_location, new_output_location)

    def make_duplicate_command(self, in_fid, out_fid, new_input_location, new_output_location):
        ## TODO: Simplify the code below
        in_chunk = new_input_location
        if(in_chunk == "stdin"):
            new_stdin = in_fid
        else:
            ## Question: Is that valid?
            new_stdin = self.stdin

        if(new_output_location == "stdout"):
            new_stdout = out_fid
        else:
            ## Question: Is that valid?
            new_stdout = self.stdout

        new_options = self.options.copy()
        if(isinstance(in_chunk, tuple)
           and len(in_chunk) == 2
           and in_chunk[0] == "option"):
            new_options[in_chunk[1]] = in_fid

        new_in_stream = [new_input_location]
        new_out_stream = [new_output_location]
        ## TODO: I probably have to do the same with output options

        new_command = Command(None, # The ast is None
                              self.command,
                              new_options,
                              new_in_stream,
                              new_out_stream,
                              self.opt_indices,
                              self.category,
                              new_stdin,
                              new_stdout)
        ## Question: Is it valid setting stdin and stdout to the stdin
        ## and stdout of the current command?
        return new_command

    ## TODO: Rename
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

class Cat(Command):
    def __init__(self, ast, command, options, in_stream, out_stream,
                 opt_indices, category, stdin=None, stdout=None):
        super().__init__(ast, command, options, in_stream, out_stream,
                         opt_indices, category, stdin=stdin, stdout=stdout)

def make_cat_node(input_file_ids, output_file_id):
    command = string_to_argument("cat")
    options = input_file_ids
    in_stream = [("option", i)  for i in range(len(input_file_ids))]
    out_stream = ["stdout"]
    stdout = output_file_id
    opt_indices = []
    category = "stateless"
    ## TODO: Fill the AST
    ast = None
    return Cat(ast, command, options, in_stream, out_stream,
               opt_indices, category, stdout=stdout)

