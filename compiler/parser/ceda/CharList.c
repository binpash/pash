#include <malloc.h>
#include <assert.h>
#include <string.h>

#include "CharList.h"


struct charList {
    struct charNode* head;
    struct charNode* last;

    int length;
};


struct charNode {
    char c;         // Translate OCaml 's : char list'

    struct charNode* next;
};


//-------------------------------------------------------------------------------------------
// char or arg_char


static void checkList (CharList myList) {
    return;

    assert (myList != NULL);

    int len = 0;
    struct charNode* cur = myList->head;

    if (myList->head == NULL) {
        assert (myList->last == NULL);
        assert (myList->length == 0);
    } else {
        while (cur != NULL) {
            if (cur->next == NULL) {
                assert (cur == myList->last);
            }

            len ++;
            cur = cur->next;
        }

        assert (len == myList->length);
    }
}


CharList newCharList (void) {
    CharList myList = malloc (sizeof (struct charList));
    assert (myList != NULL);

    myList->head = NULL;
    myList->last = NULL;
    myList->length = 0;

    checkList (myList);

    return (myList);
}


int isCharListEmpty (CharList myList) {
    assert (myList != NULL);

    if (myList->length == 0) {
        assert (myList->head == NULL);
        assert (myList->last == NULL);
    } else {
        assert (myList->head != NULL);
        assert (myList->last != NULL);
    }

    checkList (myList);

    return (myList->head == NULL);
}


int charListLength (CharList myList) {
    assert (myList != NULL);

//    checkList (myList);

    return (myList->length);
}


void charListTail (CharList myList) {
    assert (myList != NULL);

    assert (myList->length > 0);

    assert (myList->head != NULL);
    assert (myList->last != NULL);

    struct charNode* oldHead = myList->head;
    myList->head = oldHead->next;

    free (oldHead);

    if (myList->head == NULL) {
        // []
        myList->last = NULL;
    }

    myList->length = myList->length - 1;

    checkList (myList);
}


void destroyCharList (CharList myList) {
    assert (myList != NULL);

    checkList (myList);

    while (! isCharListEmpty (myList)) {
        charListTail (myList);
    }

    free (myList);
}


//-------------------------------------------------------------------------------------------
// char


void prependCharList_char (CharList myList, char newChar) {
    assert (myList != NULL);

    struct charNode* newHead = malloc (sizeof (struct charNode));
    assert (newHead != NULL);

    newHead->c = newChar;
    newHead->next = myList->head;

    if (isCharListEmpty (myList)) {
        myList->last = newHead;
    }

    myList->head = newHead;
    myList->length = myList->length + 1;
}


void appendCharList_char (CharList myList, char newChar) {
    assert (myList != NULL);

    struct charNode* newLast = malloc (sizeof (struct charNode));
    assert (newLast != NULL);

    newLast->c = newChar;
    newLast->next = NULL;

    if (isCharListEmpty (myList)) {
        myList->head = newLast;
    } else {
        myList->last->next = newLast;
    }

    myList->last = newLast;
    myList->length = myList->length + 1;
}


char charListHead_char (CharList myList) {
    assert (myList != NULL);

    assert (myList->head != NULL);
    assert (myList->last != NULL);

    return (myList->head->c);
}


char charListSecond_char (CharList myList) {
    assert (myList != NULL);

    assert (myList->head != NULL);
    assert (myList->head->next != NULL);

    return (myList->head->next->c);
}


char charListLast_char (CharList myList) {
    assert (myList != NULL);

    assert (myList->head != NULL);
    assert (myList->last != NULL);

    return (myList->last->c);
}


// Returns the first half of the list, divided by the separator (which stays in the original list).
CharList split_at (CharList myList, char separator) {
    assert (myList != NULL);

    CharList leftList = newCharList ();

    while (myList != NULL) {
        if (charListHead_char (myList) == separator) {
            break; // Ugh
        }

        appendCharList_char (leftList, charListHead_char (myList));

        charListTail (myList);
    }

    return (leftList);
}
