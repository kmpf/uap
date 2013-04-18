#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define BLOCKSIZE (4096 * 1024)

int main(int argc, char** argv)
{
    FILE* f_out = 0;
    int fd_out = 1;
    int i;
    int rsize;
    char* block = malloc(BLOCKSIZE);

    if (argc > 2)
    {
        if (strcmp(argv[argc - 2], "-o") == 0)
        {
            f_out = fopen(argv[argc - 1], "w");
            fd_out = fileno(f_out);
            argc -= 2;
        }
    }
    for (i = 1; i < argc; ++i)
    {
        FILE* f = fopen(argv[i], "r");
        while (!feof(f)) {
            rsize = fread(block, 1, BLOCKSIZE, f);
            if (rsize)
                write(fd_out, block, rsize);
        }
    }
    free(block);
    if (f_out)
        fclose(f_out);
    return 0;
}
