#include <fcntl.h>
#include <stdlib.h>
#include <stdio.h>

// cat4m reads all input files specified on the command line in chunks
// of 4M and writes them to stdout, or to another destination which can
// be specified via -o outfile (however, these two parameters must appear
// at the end of the argument list!)
//
// compile with gcc -o cat4m cat4m.c
//
// examples:
// 
// $ ./cat4m file1 file2 [...]
// $ ./cat4m file1 file2 [...] -o outfile

#define BLOCKSIZE 4194304

int main(int argc, char** argv)
{
    // write to stdout by default
    int fd_out = 1;
    int i;
    char* block = malloc(BLOCKSIZE);
    if (!block)
    {
        perror("Unable to allocate buffer");
        exit(1);
    }

    if ((argc > 2) && (strcmp(argv[argc - 2], "-o") == 0))
    {
        fd_out = open(argv[argc - 1], 
                      O_CREAT| O_WRONLY | O_TRUNC, 
                      S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH | S_IWOTH);
        if (fd_out == -1)
        {
            perror("Unable to open output file for writing");
            exit(1);
        }
        argc -= 2;
    }
    
    for (i = 1; i < argc; ++i)
    {
        int fd_in = open(argv[i], O_RDONLY);
        if (fd_in == -1)
        {
            perror("Unable to open input file for reading");
            exit(1);
        }
        while (1) {
            int rsize = read(fd_in, block, BLOCKSIZE);
            if (rsize)
            {
                int wsize = write(fd_out, block, rsize);
                if (wsize != rsize)
                {
                    perror("Error writing to output file");
                    exit(1);
                }
            }
            else
                break;
        }
    }

    if (close(fd_out) == -1)
    {
        perror("Unable to close output file");
        exit(1);
    }
    free(block);
    return 0;
}
