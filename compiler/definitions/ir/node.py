from util import *

## TODO: This might be completely obsolete now. A Node is not used in the DFG.
class Node:
    def __init__(self, ast, in_stream=[], out_stream=[],
                 category="none", stdin=None, stdout=None):
        self.ast = ast
        self.in_stream = in_stream
        self.out_stream = out_stream
        self.stdin = stdin
        self.stdout = stdout
        self.category = category

    def __repr__(self):
        output = "Node: \"{}\" in:{} out:{}".format(
            self.ast, self.stdin, self.stdout)
        return output

    ## These two commands return the flattened fileId list. Meaning
    ## that they return the children, if they exist.
    def get_input_file_ids(self):
        return [self.get_file_id(input_chunk) for input_chunk in self.in_stream]

    def get_output_file_ids(self):
        return [self.get_file_id(output_chunk) for output_chunk in self.out_stream]

    def get_number_of_inputs(self):
        return len(self.in_stream)

    def get_number_of_output(self):
        return len(self.out_stream)

    ## TODO: Rename
    def get_file_id(self, chunk):
        if (chunk == "stdout"):
            return self.stdout

        if (chunk == "stdin"):
            return self.stdin

        ## TODO: Complete this
        log(chunk)
        assert(False)

    ## TODO: Is there a way to abstract the behaviour of these two functions?
    def set_file_id(self, chunk, value):
        if (chunk == "stdout"):
            self.stdout = value
        elif (chunk == "stdin"):
            self.stdin = value
        else:
            ## TODO: Complete this
            log(chunk, value)
            assert(False)
