

## Question: Is this information adequate?
class Command:
    def __init__(self, command, stdin=None, stdout=None, options=None):
        self.command = command
        self.stdin = stdin
        self.stdout = stdout
        self.options = options
        
    def __repr__(self):
        output = "Command: {} in:{} out:{} opts:{}".format(
            self.command, self.stdin, self.stdout, self.options)
        return output

# class Arg:
#     def __init__(self, arg_char_list):
#         self.arg_char_list = arg_char_list


#     def __repr__(self):
        
