typedef struct argCharList* ArgCharList;


ArgCharList newArgCharList (void);
int isArgCharListEmpty (ArgCharList myList);
int argCharListLength (ArgCharList myList);
void* argCharListHead (ArgCharList myList);
void argCharListTail (ArgCharList myList);
void appendArgCharList (ArgCharList myList, void* new_arg_char);
void destroyArgCharList (ArgCharList myList);
