#include <malloc.h>
#include <assert.h>

#include "arg_char.h"


// Strings instead of char for convenient JSON serialization
const char* SERIALIZE_TYPE_ARG_CHAR [] = {"C", "E", "T", "A", "V", "Q", "B"};


struct arg_char_TYPE* newArgCharC (char c) {
    struct arg_char_TYPE* a = malloc (sizeof (struct arg_char_TYPE));
    assert (a != NULL);
    a->type = TYPE_ARG_CHAR_C;
    a->C.c  = c;

    return a;
}


struct arg_char_TYPE* newArgCharE (char c) {
    struct arg_char_TYPE* a = malloc (sizeof (struct arg_char_TYPE));
    assert (a != NULL);
    a->type = TYPE_ARG_CHAR_E;
    a->E.c  = c;

    return a;
}


struct arg_char_TYPE* newArgCharT (char* str) {
    struct arg_char_TYPE* a = malloc (sizeof (struct arg_char_TYPE));
    assert (a != NULL);
    a->type  = TYPE_ARG_CHAR_T;
    a->T.str = str;

    return a;
}


struct arg_char_TYPE* newArgCharA (arg_TYPE arg) {
    struct arg_char_TYPE* a = malloc (sizeof (struct arg_char_TYPE));
    assert (a != NULL);
    a->type  = TYPE_ARG_CHAR_A;
    a->A.arg = arg;

    return a;
}


struct arg_char_TYPE* newArgCharV (int var_type, int vsnul, char* str, arg_TYPE arg) {
    struct arg_char_TYPE* a = malloc (sizeof (struct arg_char_TYPE));
    assert (a != NULL);
    a->type       = TYPE_ARG_CHAR_V;
    a->V.var_type = var_type;
    a->V.vsnul    = vsnul;
    a->V.str      = str;
    a->V.arg      = arg;

    return a;
}


struct arg_char_TYPE* newArgCharQ (arg_TYPE arg) {
    struct arg_char_TYPE* a = malloc (sizeof (struct arg_char_TYPE));
    assert (a != NULL);
    a->type  = TYPE_ARG_CHAR_Q;
    a->Q.arg = arg;

    return a;
}


struct arg_char_TYPE* newArgCharB (struct t_TYPE* t) {
    struct arg_char_TYPE* a = malloc (sizeof (struct arg_char_TYPE));
    assert (a != NULL);
    a->type = TYPE_ARG_CHAR_B;
    a->B.t  = t;

    return a;
}
