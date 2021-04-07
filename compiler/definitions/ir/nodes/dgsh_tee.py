from definitions.ir.dfg_node import *

class DGSHTee(DFGNode):
    def __init__(self, inputs, outputs, com_name, com_category, com_options = [], 
                 com_redirs = [], com_assignments=[]):
        super().__init__(inputs, outputs, com_name, com_category, 
                         com_options=com_options, 
                         com_redirs=com_redirs, 
                         com_assignments=com_assignments)

def make_dgsh_tee_node(input_id, output_id):
    dgsh_tee_bin = os.path.join(config.PASH_TOP, config.config['runtime']['dgsh_tee_binary'])
    com_name = Arg(string_to_argument(dgsh_tee_bin))

    com_category = "pure"

    ## TODO: add as command line arguments
    com_options = [(2, Arg(string_to_argument("-I")))] # Eager functionality
    com_options.append((3, Arg(string_to_argument("-f")))) # use file on disk when buffer reaches maximum
    # com_options.append(4, Arg(string_to_argument("−m batch_size"))) # set the 

    return DGSHTee([input_id],
                 [output_id],
                 com_name, 
                 com_category,
                 com_options=com_options)
