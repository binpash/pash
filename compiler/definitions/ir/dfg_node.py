class DFGNode:
    ## TODO: Make inputs + outputs be structure of fids instead of list. 
    ##       This will allow us to handle comm and other commands with static inputs. 
    ##
    ## inputs : list of fids
    ## outputs : list of fids
    ## com_name : command name Arg
    ## com_options : list of arguments Arg
    ## com_redirs : list of redirections
    ## com_assignments : list of assignments
    def __init__(self, inputs, outputs, com_name, com_options = [], 
                 com_redirs = [], com_assignments=[]):
        self.inputs = inputs
        self.outputs = outputs
        self.com_name = com_name
        self.com_options = com_options
        self.com_redirs = com_redirs
        self.com_assignments = com_assignments
    

