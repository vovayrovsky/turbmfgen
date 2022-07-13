#include "random.h"

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

int GetRandi (int max, int min)
    {
    return rand() % (max - min) + min;
    }

int GetRandiMax()
    {
    return rand();
    }

double GetRandd (double max, double min)
    {
    long long buf = (long long) rand() << 31 | rand();          //This is done for more digits after decimal point
    long long max_buf = (long long) RAND_MAX << 31 | RAND_MAX;

    return min + (max - min) * (double) buf / max_buf;
    }

int InitializeRandom()
    {
    FILE* devrandom = fopen (RANDOM_DEV, "r");

    if (devrandom == NULL)
        return DEVRANDOM_OPEN_ERROR;

    int buf = 0;
    
    if (fread (&buf, sizeof (buf), 1, devrandom) != 1)
        return DEVRANDOM_READ_ERROR;

    fclose (devrandom);

    srand (buf);

    return NO_ERROR;
    }

//EOF
