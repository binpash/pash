#include "ArgCharStack.h"


#define TYPE_ARG_CHAR_C 0
#define TYPE_ARG_CHAR_E 1
#define TYPE_ARG_CHAR_T 2
#define TYPE_ARG_CHAR_A 3
#define TYPE_ARG_CHAR_V 4
#define TYPE_ARG_CHAR_Q 5
#define TYPE_ARG_CHAR_B 6


// typedef ArgCharList arg_TYPE; // arg = arg_char list
typedef ArgCharStack arg_TYPE; // arg = arg_char list


extern const char* SERIALIZE_TYPE_ARG_CHAR [];


//------------------------------------------------------------------------------------------
// arg_char = ...


//   | C of char
struct arg_char_C {
    unsigned char c;
};

//   | E of char (* escape... necessary for expansion *)
struct arg_char_E {
    unsigned char c;
};

//   | T of string option (* tilde *)
struct arg_char_T {
    char* str;
};

//   | A of arg (* arith *)
struct arg_char_A {
    arg_TYPE arg;
};

//   | V of var_type * bool (* VSNUL? *) * string * arg
struct arg_char_V {
    int var_type;
    int vsnul;
    char* str;
    arg_TYPE arg;
};

//   | Q of arg (* quoted *)
struct arg_char_Q {
    arg_TYPE arg;
};

//   | B of t (* backquote *)
struct arg_char_B {
    struct t_TYPE* t;
};

struct arg_char_TYPE {
    int type;

    union {
        struct arg_char_C C;
        struct arg_char_E E;
        struct arg_char_T T;
        struct arg_char_A A;
        struct arg_char_V V;
        struct arg_char_Q Q;
        struct arg_char_B B;
    };
};


struct arg_char_TYPE* newArgCharC (char c);
struct arg_char_TYPE* newArgCharE (char c);
struct arg_char_TYPE* newArgCharT (char* str);
struct arg_char_TYPE* newArgCharA (arg_TYPE arg);
struct arg_char_TYPE* newArgCharV (int var_type, int vsnul, char* str, arg_TYPE arg);
struct arg_char_TYPE* newArgCharQ (arg_TYPE arg);
struct arg_char_TYPE* newArgCharB (struct t_TYPE* t);
