#include <malloc.h>
#include <assert.h>
#include <string.h>

#include "CharList.h"


struct charList {
    struct charNode* head;

    int length;
};


struct charNode {
    // Form a more perfect union instead of proper polymorphism
    union {
        char c;         // Translate OCaml 's : char list'
        void* arg_char; // Translate OCaml 'arg_char'
    };

    struct charNode* next;
};


//-------------------------------------------------------------------------------------------
// char or arg_char


CharList newCharList (void) {
    CharList myList = malloc (sizeof (struct charList));
    assert (myList != NULL);

    myList->head = NULL;
    myList->length = 0;

    return (myList);
}


int isCharListEmpty (CharList myList) {
    assert (myList != NULL);

    return (myList->head == NULL);
}


int charListLength (CharList myList) {
    assert (myList != NULL);

    return (myList->length);
}


void charListTail (CharList myList) {
    assert (myList != NULL);

    assert (myList->head != NULL);

    struct charNode* head = myList->head;
    myList->head = head->next;

    free (head);

    myList->length = myList->length - 1;
}


static void prependCharList_empty (CharList myList) {
    assert (myList != NULL);

    struct charNode* newHead = malloc (sizeof (struct charNode));
    assert (newHead != NULL);

//    newHead->c = '\0';
//    newHead->arg_char = NULL;
    newHead->next = myList->head;

    myList->head = newHead;
    myList->length = myList->length + 1;
}


void destroyCharList (CharList myList) {
    assert (myList != NULL);

    while (! isCharListEmpty (myList)) {
        charListTail (myList);
    }

    free (myList);
}


//-------------------------------------------------------------------------------------------
// char


void prependCharList_char (CharList myList, char newChar) {
    assert (myList != NULL);

    prependCharList_empty (myList);

    assert (myList->head != NULL);
    myList->head->c = newChar;
}


char charListHead_char (CharList myList) {
    assert (myList != NULL);

    assert (myList->head != NULL);

    return (myList->head->c);
}


char charListSecond_char (CharList myList) {
    assert (myList != NULL);

    assert (myList->head != NULL);
    assert (myList->head->next != NULL);

    return (myList->head->next->c);
}


CharList explode (char* str) {
    assert (str != NULL);

    CharList list = newCharList ();

    for (int i = strlen (str) - 1; i >= 0; i--) {
        prependCharList_char (list, str [i]);
    }

    return (list);
}


// Kludge: reverse
/*
char* implode (CharList myList) {
    assert (myList != NULL);

    char* str = malloc (charListLength (myList) + 1);
    assert (str != NULL);

    int i = charListLength (myList);
    str [i] = '\0';

    struct charNode* head = myList->head;
    while (head != NULL) {
        i --;
        str [i] = head->c;

        head = head->next;
    }

    return (str);
}
*/


char* implode (CharList myList) {
    assert (myList != NULL);

    char* str = malloc (charListLength (myList) + 1);
    assert (str != NULL);

    int i = 0;

    struct charNode* head = myList->head;
    while (head != NULL) {
        str [i] = head->c;
        i ++;

        head = head->next;
    }

    str [i] = '\0';

    return (str);
}


/*
let rec split_at p xs =
  match xs with
  | [] -> ([],[])
  | x::xs ->
     if p x
     then ([],x::xs)
     else let (xs,ys) = split_at p xs in
          (x::xs, ys)
*/


// Returns the first half of the list, divided by the separator (which stays in the original list).
CharList split_at (CharList myList, char separator) {
    assert (myList != NULL);

    CharList leftList = newCharList ();

    while (myList != NULL) {
        if (charListHead_char (myList) == separator) {
            break; // Ugh
        }

        prependCharList_char (leftList, charListHead_char (myList));

        charListTail (myList);
    }

    return (leftList);
}


/*
int existsInCharList (CharList myList, char key) {
    assert (myList != NULL);

    struct charNode* head = myList->head;
    while (head != NULL) {
        if (head->c == key) {
            return TRUE;
        }

        head = head->next;
    }

    return FALSE;
}
*/


//-------------------------------------------------------------------------------------------
// arg_char


void prependCharList_arg_char (CharList myList, void* new_arg_char) {
    assert (myList != NULL);

    prependCharList_empty (myList);

    assert (myList->head != NULL);
    myList->head->arg_char = new_arg_char;
}


void* charListHead_arg_char (CharList myList) {
    assert (myList != NULL);

    assert (myList->head != NULL);

    return (myList->head->arg_char);
}
