#include <malloc.h>
#include <assert.h>
#include <string.h>

#include "ArgCharList.h"


struct argCharList {
    struct argCharNode* head;
    struct argCharNode* last;

    int length;
};


struct argCharNode {
    void* arg_char; // Translate OCaml 'arg_char'

    struct argCharNode* next;
};


ArgCharList newArgCharList (void) {
    ArgCharList myList = malloc (sizeof (struct argCharList));
    assert (myList != NULL);

    myList->head = NULL;
    myList->last = NULL;
    myList->length = 0;

    return (myList);
}


int isArgCharListEmpty (ArgCharList myList) {
    assert (myList != NULL);

    if (myList->length == 0) {
        assert (myList->head == NULL);
        assert (myList->last == NULL);
    } else {
        assert (myList->head != NULL);
        assert (myList->last != NULL);
    }

    return (myList->head == NULL);
}


int argCharListLength (ArgCharList myList) {
    assert (myList != NULL);

    return (myList->length);
}


void* argCharListHead (ArgCharList myList) {
    assert (myList != NULL);

    assert (myList->head != NULL);

    return (myList->head->arg_char);
}


void argCharListTail (ArgCharList myList) {
    assert (myList != NULL);

    assert (myList->length > 0);

    assert (myList->head != NULL);

    struct argCharNode* oldHead = myList->head;
    myList->head = oldHead->next;

    free (oldHead);

    if (myList->head == NULL) {
        myList->last = NULL;
    } else if (myList->head->next == NULL) {
        myList->last = myList->head;
    }

    myList->length = myList->length - 1;
}


// Add to end of the list.
void appendArgCharList (ArgCharList myList, void* new_arg_char) {
    assert (myList != NULL);

    struct argCharNode* newLast = malloc (sizeof (struct argCharNode));
    assert (newLast != NULL);

    newLast->arg_char = new_arg_char;
    newLast->next = NULL;

    if (myList->head == NULL) {
        myList->head = newLast;
        myList->last = newLast;
    } else {
        myList->last->next = newLast;
        myList->last = newLast;
    }

    myList->length = myList->length + 1;
}


void destroyArgCharList (ArgCharList myList) {
    assert (myList != NULL);

    while (! isArgCharListEmpty (myList)) {
        argCharListTail (myList);
    }

    free (myList);
}
