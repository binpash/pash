typedef struct charList* CharList;


CharList newCharList (void);
int isCharListEmpty (CharList myList);
int charListLength (CharList myList);
void charListTail (CharList myList);
void charListInit (CharList myList);
void destroyCharList (CharList myList);

void prependCharList_char (CharList myList, char newChar);
void appendCharList_char (CharList myList, char newChar);
char charListHead_char (CharList myList);
char charListSecond_char (CharList myList);
CharList split_at (CharList myList, char separator);
// int existsInCharList (CharList myList, char key);
