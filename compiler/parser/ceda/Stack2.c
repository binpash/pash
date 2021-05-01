#include <string.h>
#include <malloc.h>
#include <assert.h>

#include "Stack.h"


#define TRUE  1
#define FALSE 0


#define INIT_STACK_SIZE 32


struct stack {
    unsigned int allocSize;
    unsigned int usedSize;

    char* items;
};


Stack newStack (void) {
    Stack myStack = malloc (sizeof (struct stack));
    assert (myStack != NULL);

    myStack->allocSize = INIT_STACK_SIZE;
    myStack->usedSize  = 0;

    myStack->items = malloc (sizeof (char) * myStack->allocSize);
    assert (myStack->items != NULL);

    return myStack;
}



int isStackEmpty (Stack myStack) {
    assert (myStack != NULL);

    return (myStack->usedSize == 0);
}


void pushStack (Stack myStack, char item) {
    assert (myStack != NULL);

    if (myStack->usedSize >= myStack->allocSize) {
        myStack->allocSize = myStack->allocSize * 2;
        myStack->items = realloc (myStack->items, sizeof (char) * myStack->allocSize);
        assert (myStack->items != NULL);
    }

    myStack->items [myStack->usedSize] = item;
    myStack->usedSize ++;
}


char popStack (Stack myStack) {
    assert (myStack != NULL);

    assert (myStack->usedSize > 0);

    myStack->usedSize --;
    char item = myStack->items [myStack->usedSize];
    // Don't bother shrinking stack

    return item;
}


char topStack (Stack myStack) {
    assert (myStack != NULL);

    assert (myStack->usedSize > 0);

    char item = myStack->items [myStack->usedSize - 1];

    return (item);
}


char secondTopStack (Stack myStack) {
    assert (myStack != NULL);

    assert (myStack->usedSize > 1);

    char item = myStack->items [myStack->usedSize - 2];

    return (item);
}


unsigned int getStackSize (Stack myStack) {
    assert (myStack != NULL);

    return (myStack->usedSize);
}


void destroyStack (Stack myStack) {
    assert (myStack != NULL);

    assert (myStack->items != NULL);
    free (myStack->items);

    free (myStack);
}


char* serializeStack (Stack myStack) {
    assert (myStack != NULL);

    char* str = malloc (myStack->usedSize + 1);
    assert (str != NULL);

    memcpy (str, myStack->items, sizeof (char) * myStack->usedSize);
    str [myStack->usedSize] = '\0';

    return str;
}


int existsInStack (Stack myStack, char key) {
    assert (myStack != NULL);

    for (int i = 0; i < myStack->usedSize; i++) {
        if (myStack->items [i] == key) {
            return TRUE;
        }
    }

    return FALSE;
}


// Generic implementation
// Deliberately reversed
Stack explode_rev0 (char* str) {
    assert (str != NULL);

    Stack list = newStack ();

    int len = strlen (str);
    for (int i = len - 1; i >= 0; i--) {
//    for (int i = 0; i < len; i++) {
        pushStack (list, str [i]);
    }

    return (list);
}


// Deliberately reversed
Stack explode_rev (char* str) {
    assert (str != NULL);

    Stack myStack = newStack ();
    assert (myStack != NULL);

    int len = strlen (str);

    if (myStack->allocSize < len) {
        myStack->allocSize = len;
        myStack->items = realloc (myStack->items, sizeof (char) * myStack->allocSize);
        assert (myStack->items != NULL);
    }

    for (int i = 0; i < len; i++) {
        myStack->items [len - 1 - i] = str [i];
    }

    myStack->usedSize = len;

    return (myStack);
}


char* implode0 (Stack myList) {
    assert (myList != NULL);

    char* str = malloc (getStackSize (myList) + 1);
    assert (str != NULL);

    int i = 0;
    while (getStackSize (myList) > 0) {
        str [i] = popStack (myList);

        i ++;
    }

    str [i] = '\0';

    destroyStack (myList);

    return (str);
}


char* implode_rev (Stack myStack) {
    assert (myStack != NULL);

    int len = getStackSize (myStack);
    char* str = malloc (sizeof (char) * (len + 1));
    assert (str != NULL);

    for (int i = 0; i < len; i++) {
        str [i] = myStack->items [len - 1 - i];
    }
    str [len] = '\0';

    return (str);
}


// Returns the first half of the list, divided by the separator (which stays in the original list).
Stack split_at (Stack myStack, char separator) {
    assert (myStack != NULL);

    Stack leftStack = newStack ();

    while (getStackSize (myStack) > 0) {
        if (topStack (myStack) == separator) {
            break; // Ugh
        }

        pushStack (leftStack, popStack (myStack));
    }

    Stack revLeftStack = newStack ();
    while (getStackSize (leftStack) > 0) {
        pushStack (revLeftStack, popStack (leftStack));
    }

    return (revLeftStack);
}


/*



// Returns the first half of the list, divided by the separator (which stays in the original list).
Stack split_at (Stack myStack, char separator) {
    assert (myStack != NULL);

    Stack leftStack = newStack ();

    while (getStackSize (myStack) > 0) {
        if (topStack (myStack) == separator) {
            break; // Ugh
        }

        pushStack (leftStack, popStack (myStack));
    }

    return (leftStack);
}
*/
