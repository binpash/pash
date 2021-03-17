typedef struct stack* Stack;

Stack newStack (void);
int isStackEmpty (Stack myStack);
void pushStack (Stack myStack, char newChar);
char popStack (Stack myStack);
char topStack (Stack myStack);
char* serializeStack (Stack myStack);
int existsInStack (Stack myStack, char key);
void destroyStack (Stack myStack);
