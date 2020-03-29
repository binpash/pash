import parse
import json
import wrapper
import config

print("Parse binary {}".format(config.PARSER_BINARY))

ast = parse.shell_file_to_ir("./t1.sh")
wrapper.rewrite_ast(ast)
print("Wrapped {} commands".format(wrapper.get_results()[0]))
parse.ir_to_shell(ast, "./t2.sh")
