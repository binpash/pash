#ifndef AST2_H
#define AST2_H


#include "ArgCharStack.h"


#define TRUE  1
#define FALSE 0


//------------------------------------------------------------------------------------------


#define STACK_CTLVar 100
#define STACK_CTLAri 101
#define STACK_CTLQuo 102


//------------------------------------------------------------------------------------------


#define TYPE_T_COMMAND    0
#define TYPE_T_PIPE       1
#define TYPE_T_REDIR      2
#define TYPE_T_BACKGROUND 3
#define TYPE_T_SUBSHELL   4
#define TYPE_T_AND        5
#define TYPE_T_OR         6
#define TYPE_T_NOT        7
#define TYPE_T_SEMI       8
#define TYPE_T_IF         9
#define TYPE_T_WHILE      10
#define TYPE_T_FOR        11
#define TYPE_T_CASE       12
#define TYPE_T_DEFUN      13

extern const char* SERIALIZE_TYPE_T [];


/*
    var_type =
       | Normal
       | Minus
       | Plus
       | Question
       | Assign
       | TrimR
       | TrimRMax
       | TrimL
       | TrimLMax
       | Length
*/
#define VAR_TYPE_NORMAL    0x0
#define VAR_TYPE_MINUS     0x2
#define VAR_TYPE_PLUS      0x3
#define VAR_TYPE_QUESTION  0x4
#define VAR_TYPE_ASSIGN    0x5
#define VAR_TYPE_TRIMR     0x6
#define VAR_TYPE_TRIMRMAX  0x7
#define VAR_TYPE_TRIML     0x8
#define VAR_TYPE_TRIMLMAX  0x9
#define VAR_TYPE_LENGTH    0xA

extern const char* SERIALIZE_VAR_TYPE [];


/*
   redirection =
   | File of redir_type * int * arg
   | Dup of dup_type * int * arg
   | Heredoc of heredoc_type * int * arg
*/
// Don't mix up with REDIR_TYPE_*!
#define REDIRECTION_TYPE_FILE    0x0
#define REDIRECTION_TYPE_DUP     0x1
#define REDIRECTION_TYPE_HEREDOC 0x2

extern const char* SERIALIZE_REDIRECTION_TYPE [];


// redir_type = To | Clobber | From | FromTo | Append
#define REDIR_TYPE_TO      0x0
#define REDIR_TYPE_CLOBBER 0x1
#define REDIR_TYPE_FROM    0x2
#define REDIR_TYPE_FROMTO  0x3
#define REDIR_TYPE_APPEND  0x4

extern const char* SERIALIZE_REDIR_TYPE [];


// dup_type = ToFD | FromFD
#define DUP_TYPE_TOFD   0x0
#define DUP_TYPE_FROMFD 0x1

extern const char* SERIALIZE_DUP_TYPE [];


// heredoc_type = Here | XHere (* for when in a quote... not sure when this comes up *)
#define HEREDOC_TYPE_HERE  0x0
#define HEREDOC_TYPE_XHERE 0x1

extern const char* SERIALIZE_HEREDOC_TYPE [];


// Duplicates arg_char.h
typedef ArgCharStack arg_TYPE; // arg = arg_char list


//------------------------------------------------------------------------------------------
// SIMPLE LIST TYPES


struct t_list {
    struct t_TYPE* t;
    struct t_list* next;
};

struct assign_list {
    struct assign_TYPE* assign;
    struct assign_list* next;
};

struct redirectionList {
    struct redirection_TYPE* redir;
    struct redirectionList* next;
};

// args = arg list
// Note that 'arg_TYPE' is typedef'ed above as a CharList
struct args_TYPE {
    arg_TYPE arg;
    struct args_TYPE* next;
};

struct case_list {
    struct case_TYPE* casey;
    struct case_list* next;
};


//------------------------------------------------------------------------------------------


// assign = string * arg
struct assign_TYPE {
    char* string;
    arg_TYPE arg;
};


// | File of redir_type * int * arg
//
// Not to be mistaken with <stdio.h>'s FILE
struct file_TYPE {
    int redir_type;
    int fd;
    arg_TYPE a;
};


// | Dup of dup_type * int * arg
struct dup_TYPE {
    int dup_type;
    int fd;
    arg_TYPE tgt;
};


// | Heredoc of heredoc_type * int * arg
struct heredoc_TYPE {
    int heredoc_type;
    int fd;
    arg_TYPE a;
};


/*
   redirection =
   | File of redir_type * int * arg
   | Dup of dup_type * int * arg
   | Heredoc of heredoc_type * int * arg
*/
struct redirection_TYPE {
    int type;

    union {
        struct file_TYPE file;
        struct dup_TYPE dup;
        struct heredoc_TYPE heredoc;
    };
};


// case = { cpattern : arg list; cbody : t }
//
// Note: the hash table is useful only for JSON serialization/deserialization
// purposes; we don't store it that way internally.
struct case_TYPE {
    struct args_TYPE* cpattern;
    struct t_TYPE* cbody;
};


//------------------------------------------------------------------------------------------
// type t = ...


//  | Command of linno * assign list * args * redirection list (* assign, args, redir *)
struct Command_TYPE {
    unsigned int linno;
    struct assign_list* assign;
    struct args_TYPE* args;
    struct redirectionList* redirect;
};


//  | Pipe of bool * t list (* background?, commands *)
struct Pipe_TYPE {
    int background;
    struct t_list* spill;
};


//  | Redir of linno * t * redirection list
struct Redir_TYPE {
    unsigned int linno;
    struct t_TYPE* t;
    struct redirectionList* redirect;
};


//  | Background of linno * t * redirection list 
struct Background_TYPE {
    unsigned int linno;
    struct t_TYPE* t;
    struct redirectionList* redirect;
};


//  | Subshell of linno * t * redirection list
struct Subshell_TYPE {
    unsigned int linno;
    struct t_TYPE* t;
    struct redirectionList* redirect;
};


//  | And of t * t
struct And_TYPE {
    struct t_TYPE* left;
    struct t_TYPE* right;
};


//  | Or of t * t
struct Or_TYPE {
    struct t_TYPE* left;
    struct t_TYPE* right;
};


//  | Not of t
struct Not_TYPE {
    struct t_TYPE* t;
};


//  | Semi of t * t
struct Semi_TYPE {
    struct t_TYPE* left;
    struct t_TYPE* right;
};


//  | If of t * t * t (* cond, then, else *)
struct If_TYPE {
    struct t_TYPE* test;
    struct t_TYPE* ifpart;
    struct t_TYPE* elsepart;
};


//  | While of t * t (* test, body *) (* until encoded as a While . Not *)
struct While_TYPE {
    struct t_TYPE* test;
    struct t_TYPE* body;
};


//  | For of linno * arg * t * string (* args, body, var *)
struct For_TYPE {
    unsigned int linno;
    arg_TYPE arg;
    struct t_TYPE* body;
    char* var;
};


// | Case of linno * arg * case list
struct Case_TYPE {
    unsigned int linno;
    arg_TYPE arg;
    struct case_list* cases;
};


//  | Defun of linno * string * t (* name, body *)
struct Defun_TYPE {
    unsigned int linno;
    char* name;
    struct t_TYPE* body;
};


struct t_TYPE {
    int type;

    union {
        struct Command_TYPE Command;
        struct Pipe_TYPE Pipe;
        struct Redir_TYPE Redir;
        struct Background_TYPE Background;
        struct Subshell_TYPE Subshell;
        struct And_TYPE And;
        struct Or_TYPE Or;
        struct Not_TYPE Not;
        struct Semi_TYPE Semi;
        struct If_TYPE If;
        struct While_TYPE While;
        struct For_TYPE For;
        struct Case_TYPE Case;
        struct Defun_TYPE Defun;
    };
};


//------------------------------------------------------------------------------------------


struct t_TYPE* of_node (union node* n);


#endif
