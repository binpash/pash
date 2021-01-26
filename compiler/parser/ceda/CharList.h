typedef struct charList* CharList;


CharList newCharList (void);
int isCharListEmpty (CharList myList);
int charListLength (CharList myList);
void charListTail (CharList myList);
void destroyCharList (CharList myList);

void prependCharList_char (CharList myList, char newChar);
char charListHead_char (CharList myList);
char charListSecond_char (CharList myList);
CharList explode (char* str);
char* implode (CharList myList);
CharList split_at (CharList myList, char separator);
// int existsInCharList (CharList myList, char key);

void prependCharList_arg_char (CharList myList, void* new_arg_char);
void* charListHead_arg_char (CharList myList);
