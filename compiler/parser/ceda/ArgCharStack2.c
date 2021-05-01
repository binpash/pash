#include <string.h>
#include <malloc.h>
#include <assert.h>

#include "ArgCharStack.h"


#define TRUE  1
#define FALSE 0


#define INIT_STACK_SIZE 32


struct argCharStack {
    unsigned int allocSize;
    unsigned int usedSize;

    void** items;
};


ArgCharStack newArgCharStack (void) {
    ArgCharStack myArgCharStack = malloc (sizeof (struct argCharStack));
    assert (myArgCharStack != NULL);

    myArgCharStack->allocSize = INIT_STACK_SIZE;
    myArgCharStack->usedSize  = 0;

    myArgCharStack->items = malloc (sizeof (void*) * myArgCharStack->allocSize);
    assert (myArgCharStack->items != NULL);

    return myArgCharStack;
}



int isArgCharStackEmpty (ArgCharStack myArgCharStack) {
    assert (myArgCharStack != NULL);

    return (myArgCharStack->usedSize == 0);
}


void pushArgCharStack (ArgCharStack myArgCharStack, void* item) {
    assert (myArgCharStack != NULL);

    if (myArgCharStack->usedSize >= myArgCharStack->allocSize) {
        myArgCharStack->allocSize = myArgCharStack->allocSize * 2;
        myArgCharStack->items = realloc (myArgCharStack->items, sizeof (void*) * myArgCharStack->allocSize);
        assert (myArgCharStack->items != NULL);
    }

    myArgCharStack->items [myArgCharStack->usedSize] = item;
    myArgCharStack->usedSize ++;
}


void* popArgCharStack (ArgCharStack myArgCharStack) {
    assert (myArgCharStack != NULL);

    assert (myArgCharStack->usedSize > 0);

    myArgCharStack->usedSize --;
    void* item = myArgCharStack->items [myArgCharStack->usedSize];
    // Don't bother shrinking stack

    return item;
}


void* topArgCharStack (ArgCharStack myArgCharStack) {
    assert (myArgCharStack != NULL);

    assert (myArgCharStack->usedSize > 0);

    void* item = myArgCharStack->items [myArgCharStack->usedSize - 1];

    return (item);
}


void* secondTopArgCharStack (ArgCharStack myArgCharStack) {
    assert (myArgCharStack != NULL);

    assert (myArgCharStack->usedSize > 1);

    void* item = myArgCharStack->items [myArgCharStack->usedSize - 2];

    return (item);
}


unsigned int getArgCharStackSize (ArgCharStack myArgCharStack) {
    assert (myArgCharStack != NULL);

    return (myArgCharStack->usedSize);
}


void reverseArgCharStack (ArgCharStack myArgCharStack) {
    assert (myArgCharStack != NULL);

    assert (myArgCharStack->items != NULL);

    for (int i = 0; i < myArgCharStack->usedSize / 2; i++) {
        void* left = myArgCharStack->items [i];

//        fprintf (stderr, "%d %p, %d %p\n", i, myArgCharStack->items [i],
//                 myArgCharStack->usedSize - 1 - i, myArgCharStack->items [myArgCharStack->usedSize - 1 - i]);
        myArgCharStack->items [i] = myArgCharStack->items [myArgCharStack->usedSize - 1 - i];

        myArgCharStack->items [myArgCharStack->usedSize - 1 - i] = left;
    }
}


void destroyArgCharStack (ArgCharStack myArgCharStack) {
    assert (myArgCharStack != NULL);

    assert (myArgCharStack->items != NULL);
    free (myArgCharStack->items);

    free (myArgCharStack);
}
