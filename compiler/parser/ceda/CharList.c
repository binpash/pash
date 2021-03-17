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
    // Form a more perfect union instead of proper polymorphism
    union {
        char c;         // Translate OCaml 's : char list'
        void* arg_char; // Translate OCaml 'arg_char'
    };

    struct charNode* next;
    struct charNode* prev;
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
    } else {
        // head : rest
        myList->head->prev = NULL;
    }

    myList->length = myList->length - 1;

    checkList (myList);
}


// Haskell-style: "Return all the elements of a list except the last one."
void charListInit (CharList myList) {
    assert (myList != NULL);

    assert (myList->head != NULL);
    assert (myList->last != NULL);

    struct charNode* oldLast = myList->last;
    myList->last = oldLast->prev;

    free (oldLast);

    if (myList->last == NULL) {
        myList->head = NULL;
    } else {
        myList->last->next = NULL;
    }

    myList->length = myList->length - 1;

    checkList (myList);
}


// Add to the start of the list
static void prependCharList_empty (CharList myList) {
    assert (myList != NULL);

    struct charNode* newHead = malloc (sizeof (struct charNode));
    assert (newHead != NULL);

    newHead->prev = NULL;
    newHead->next = myList->head;

    if (isCharListEmpty (myList)) {
        myList->last = newHead;
    } else {
        myList->head->prev = newHead;
    }

    myList->head = newHead;
    myList->length = myList->length + 1;

    checkList (myList);
}


// Add to the end of the list
static void appendCharList_empty (CharList myList) {
    assert (myList != NULL);

    struct charNode* newLast = malloc (sizeof (struct charNode));
    assert (newLast != NULL);

    newLast->next = NULL;

    if (isCharListEmpty (myList)) {
        newLast->prev = NULL;
        myList->head = newLast;
    } else {
        newLast->prev = myList->last;
        myList->last->next = newLast;
    }

    myList->last = newLast;
    myList->length = myList->length + 1;

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

    prependCharList_empty (myList);

    assert (myList->head != NULL);
    assert (myList->last != NULL);
    myList->head->c = newChar;
}


void appendCharList_char (CharList myList, char newChar) {
    assert (myList != NULL);

    appendCharList_empty (myList);

    assert (myList->head != NULL);
    assert (myList->last != NULL);
    myList->last->c = newChar;
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


/*
CharList explode (char* str) {
    assert (str != NULL);

    CharList list = newCharList ();

    for (int i = strlen (str) - 1; i >= 0; i--) {
        prependCharList_char (list, str [i]);
    }

    return (list);
}
*/


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

/*
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
*/

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

        appendCharList_char (leftList, charListHead_char (myList));

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
    assert (myList->last != NULL);
    myList->head->arg_char = new_arg_char;
}


void appendCharList_arg_char (CharList myList, void* new_arg_char) {
    assert (myList != NULL);

    appendCharList_empty (myList);

    assert (myList->head != NULL);
    assert (myList->last != NULL);
    myList->last->arg_char = new_arg_char;
}


void* charListHead_arg_char (CharList myList) {
    assert (myList != NULL);

    assert (myList->head != NULL);
    assert (myList->last != NULL);

    return (myList->head->arg_char);
}


void* charListLast_arg_char (CharList myList) {
    assert (myList != NULL);

    assert (myList->head != NULL);
    assert (myList->last != NULL);

    return (myList->last->arg_char);
}
