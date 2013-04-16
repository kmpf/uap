#include <stdio.h>
#include <stdlib.h>

#define BLOCKSIZE (4096 * 1024)

int main(int argc, char** argv)
{
    int i;
    int rsize;
    char* block = malloc(BLOCKSIZE);
    for (i = 1; i < argc; ++i)
    {
        FILE* f = fopen(argv[i], "r");
        while (!feof(f)) {
            rsize = fread(block, 1, BLOCKSIZE, f);
            if (rsize)
                write(1, block, rsize);
        }
    }
    free(block);
    return 0;
}
