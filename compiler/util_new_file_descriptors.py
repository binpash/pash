from util import log
from definitions.ir.file_id import FileId
from definitions.ir.resource import FileResource, Resource
from definitions.ir.arg import Arg
# imports from annotation framework
import sys
sys.path.insert(1, "/home/felix/git-repos/MIT/annotations")
# for typing
# --
# for use
from datatypes_new.BasicDatatypes import FileName, StdDescriptor


def resource_from_file_descriptor(file_descriptor) -> Resource:
    if isinstance(file_descriptor, FileName):
        arg = Arg.string_to_arg(file_descriptor.get_name())
        log(f'filedes name: {file_descriptor.get_name()}')
        log(f'filedes name type: {type(file_descriptor.get_name())}')
        log(f'arg: {arg}')
        return FileResource(Arg.string_to_arg(file_descriptor.get_name()))
    elif isinstance(file_descriptor, StdDescriptor):
        assert(False)
        #     ## TODO: Make this be a subtype of Resource
        #     if (opt_or_fd == "stdin"):
        #         resource = ("fd", 0)
        #     elif (opt_or_fd == "stdout"):
        #         resource = ("fd", 1)
        #     elif (opt_or_fd == "stderr"):
        #         resource = ("fd", 2)
        #     else:
        #         raise NotImplementedError()
        #     resource = FileDescriptorResource(resource)
    else:
        assert(False)
    # unreachable
