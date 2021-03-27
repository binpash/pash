#include <malloc.h>
#include <assert.h>

#include "Stack.h"


#define TRUE  1
#define FALSE 0


struct node {
    char c;
    struct node* next;
};


struct stack {
    int size;
    struct node* head;
};


Stack newStack (void) {
    Stack myStack = malloc (sizeof (struct stack));
    assert (myStack != NULL);

    myStack->size = 0;
    myStack->head = NULL;

    return myStack;
}



int isStackEmpty (Stack myStack) {
    assert (myStack != NULL);

    return (myStack->size == 0);
}


void pushStack (Stack myStack, char newChar) {
    assert (myStack != NULL);

    myStack->size ++;

    struct node* myNode = malloc (sizeof (struct node));
    assert (myNode != NULL);

    myNode->c = newChar;
    myNode->next = myStack->head;

    myStack->head = myNode;
}


char popStack (Stack myStack) {
    assert (myStack != NULL);

    myStack->size --;
    assert (myStack->head != NULL);

    struct node* oldHead = myStack->head;

    char myChar = myStack->head->c;
    myStack->head = myStack->head->next;

    free (oldHead);

    return myChar;
}



char topStack (Stack myStack) {
    assert (myStack != NULL);

    assert (myStack->head != NULL);

    return (myStack->head->c);
}


void destroyStack (Stack myStack) {
    assert (myStack != NULL);

    while (! isStackEmpty (myStack)) {
        popStack (myStack);
    }

    free (myStack);
}


char* serializeStack (Stack myStack) {
    assert (myStack != NULL);

    char* str = malloc (myStack->size + 1);
    assert (str != NULL);

    int i = 0;
    struct node* cur = myStack->head;
    while (cur != NULL) {
        str [i] = cur->c;

        i ++;
        cur = cur->next;
    }
    str [i] = '\0';

    return str;
}


int existsInStack (Stack myStack, char key) {
    assert (myStack != NULL);

    struct node* head = myStack->head;
    while (head != NULL) {
        if (head->c == key) {
            return TRUE;
        }

        head = head->next;
    }

    return FALSE;
}
