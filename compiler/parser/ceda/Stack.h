typedef struct stack* Stack;

Stack newStack (void);
int isStackEmpty (Stack myStack);
void pushStack (Stack myStack, char newChar);
char popStack (Stack myStack);
char topStack (Stack myStack);
char secondTopStack (Stack myStack);
unsigned int getStackSize (Stack myStack);
void destroyStack (Stack myStack);

char* serializeStack (Stack myStack);
int existsInStack (Stack myStack, char key);

Stack explode_rev (char* str);
char* implode_rev (Stack myStack);

Stack split_at (Stack myStack, char separator);
