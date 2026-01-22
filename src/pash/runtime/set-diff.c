#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int comp (const void * elem1, const void * elem2) 
{
    char f = *((char*)elem1);
    char s = *((char*)elem2);
    if (f > s) return  1;
    if (f < s) return -1;
    return 0;
}

int main(int argc, char* argv[]) {
    // arg#1 -> old set
    // arg#2 -> new set
    if (argc != 3) {
        printf("ERROR: wrong input!\n");
        exit(1);
    }

    // Returns the - and + difference of the two sets of sorted characters.
    //
    // First ensures that characters are already sorted
    //
    // It requires that there are no duplicates
    int fromLen = strlen(argv[1]);
    int toLen = strlen(argv[2]);

    qsort(argv[1], fromLen, sizeof(char), comp);
    qsort(argv[2], toLen, sizeof(char), comp);

    int maxLen = (fromLen < toLen) ? toLen : fromLen;

    char* toRemove = malloc((maxLen + 1) * sizeof(char));
    char* toAdd = malloc((maxLen + 1) * sizeof(char));

    int iFrom = 0; 
    int iTo = 0;
    int iRemove = 0; 
    int iAdd = 0;

    while (iFrom < fromLen && iTo < toLen )
    {
        // We do not add or remove -s and -c
        if (argv[1][iFrom] == 'c' || argv[1][iFrom] == 's')
        {
            iFrom++;
        }
        else if (argv[2][iTo] == 'c' || argv[2][iTo] == 's')
        {
            iTo++;
        }
        else if (argv[1][iFrom] < argv[2][iTo])
        {
            // Only exists in iFrom, therefore we need to remove
            toRemove[iRemove] = argv[1][iFrom];
            iRemove ++;
            iFrom ++;
        }
        else if (argv[1][iFrom] == argv[2][iTo])
        {
            // Don't add or remove anything!
            iFrom ++;
            iTo ++;
        }
        else
        {
            // Only exists in iTo, therefore we need to add
            toAdd[iAdd] = argv[2][iTo];
            iAdd ++;
            iTo ++;
        }
    }

    // If there is something left in iFrom
    while (iFrom < fromLen)
    {
        if (argv[1][iFrom] != 'c' && argv[1][iFrom] != 's')
        {
            toRemove[iRemove] = argv[1][iFrom];
            iRemove ++;
        }
        iFrom ++;
    }

    // If there is something left in iTo
    while (iTo < toLen)
    {
        if (argv[2][iTo] != 'c' && argv[2][iTo] != 's')
        {
            toAdd[iAdd] = argv[2][iTo];
            iAdd ++;
        }
        iTo ++;
    }

    toRemove[iRemove] = '\0';
    toAdd[iAdd] = '\0';

    printf("%s,%s\n", toRemove, toAdd);

    return 0;
}
