typedef struct argCharStack* ArgCharStack;

ArgCharStack newArgCharStack (void);
int isArgCharStackEmpty (ArgCharStack myArgCharStack);
void pushArgCharStack (ArgCharStack myArgCharStack, void* newChar);
void* popArgCharStack (ArgCharStack myArgCharStack);
void* topArgCharStack (ArgCharStack myArgCharStack);
void* secondTopArgCharStack (ArgCharStack myArgCharStack);
unsigned int getArgCharStackSize (ArgCharStack myArgCharStack);
void reverseArgCharStack (ArgCharStack myArgCharStack);
void destroyArgCharStack (ArgCharStack myArgCharStack);
