from ctypes import *

# ../libdash/src/memalloc.h:struct
# struct stackmark {
#     struct stack_block *stackp;
#     char *stacknxt;
#     size_t stacknleft;
# };

class STACKMARK (Structure):
    _fields_ = [("stack_block", c_void_p),
                ("stacknxt",    c_char_p),
                ("stacknleft",  c_size_t)];

cdll.LoadLibrary ("../libdash/src/.libs/libdash.so");
