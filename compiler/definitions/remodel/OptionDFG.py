from definitions.remodel import IOVar

class OptionDFG:

    def __init__(self, name: str, option_arg_var: IOVar) -> None:
        self.option_name = name
        self.option_arg_var: IOVar = option_arg_var

    def get_name(self) -> str:
        return self.option_name
