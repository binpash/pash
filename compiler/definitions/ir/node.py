from util import *
from union_find import *

## Question: Is this information adequate?
##
## TODO: What other information should a node of the IR contain?
## (other redirections possibly?).
##
## (LATER) TODO: Replace all the file references in IR nodes with new
## Identifiers that we make. IN order to do this, one has to be able
## to find these file arguments (based on the analysis that we will
## do).
##
## A node represents an abstract program that our system can
## distribute. At the moment, that is a program with one input and one
## output stream. Input and output streams are shown as a list of
## either options or standard channels (such as stdin, stdout,
## stderr).
##
## Nodes also have a category, which shows whether they can be
## parallelized on their input stream or not.
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
        print(chunk)
        assert(False)

    ## TODO: Is there a way to abstract the behaviour of these two functions?
    def set_file_id(self, chunk, value):
        if (chunk == "stdout"):
            self.stdout = value
        elif (chunk == "stdin"):
            self.stdin = value
        else:
            ## TODO: Complete this
            print(chunk, value)
            assert(False)

    def find_file_id_in_in_stream(self, fileId):
        return self.find_file_id_in_stream(fileId, self.in_stream)

    def find_file_id_in_out_stream(self, fileId):
        return self.find_file_id_in_stream(fileId, self.out_stream)

    def find_file_id_in_stream(self, file_id, stream):
        index = 0
        for chunk in stream:
            chunk_file_id = Find(self.get_file_id(chunk))
            if(Find(file_id) == chunk_file_id):
                return index
            index += 1
        return None
