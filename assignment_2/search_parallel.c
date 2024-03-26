#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <stdbool.h>
#include "mpi.h"

// #define ROOT 0
#define base 256
#define prime 101

/** search()
 *  This returns the index of every match that the pattern is found
 *  in the text. This search algorithmn uses the parallelized version
 *  of the Rabin-Karp string searching algorithm using a rolling hash
 *  To find the part that each process gets we can find the amount 
 *  of indexes, where t = text length, p = pattern length
 *  and proc = number of processors
 *  n = (t-(p-1)) / proc
 *  we only need to check a maximum of t-(p-1) characters as the 
 *  pattern will be less than the charcters remaining after that point.
 *  As a result, each processor only needs to check n indexes
 *  where the indexes for each processor depends on its rank R
 *  where [R*n, (R+1)*n -1] and the last rank checks [R*n, t] indexes
 *  however, we need to pass [R*n, ((R+1)*n-1) + (p-1)] to be able to
 *  check the indexes from [((R+1)*n-1)-p, (R+1)*n-1)]
*/
void search(char* pattern, char* text, int pattern_length, int pattern_value, int hash_value);

int main(int argc, char *argv[]) {
    struct stat statbuf;
    int id;                 /* rank of executing process    */
    int p;                  /* number of processes          */
    char* pattern;
    char* file_name;
    int file_size;
    char text[file_size];
    FILE *file;
    int pattern_value;
    int hash_value;
    MPI_Status status ; /* result of MPI_Recv */
    char *elements;                 /* elements to send*/
    int pattern_length = strlen(pattern);
    int n = (file_size-(pattern_length-1)) / p;

    MPI_Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &id );
    MPI_Comm_size (MPI_COMM_WORLD, &p);

    const int ROOT = p - 1;

    MPI_Bcast(&pattern_value, 1, MPI_INT, ROOT, MPI_COMM_WORLD); 
    MPI_Bcast(&hash_value, 1, MPI_INT, ROOT, MPI_COMM_WORLD); 
    MPI_Bcast(&file_size, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
    MPI_Barrier (MPI_COMM_WORLD);
    
    /* Check usage */
    if (argc != 3) {
        if(ROOT == id){
            printf("Usage: %s <pattern> <file>\n", argv[0]);
        }
        MPI_Finalize(); 
        exit(1);
    }

    if (ROOT == id) {
        /* Call stat () to fill statbuf struct */ 
        if (lstat(argv[2], &statbuf ) == -1 ) {
            if(ROOT == id){
                printf("Could not stat file %s \n" , argv[2]);
            }
            MPI_Finalize(); 
            exit(1);
        }

        pattern = argv[1];
        file_name = argv[2];
        file_size = statbuf.st_size;

        file = fopen(file_name, "r");
        if (file == NULL) {
            if(ROOT == id){
                printf("Error opening file");
            }
            MPI_Finalize(); 
            exit(1);
        }

        /**
         * send the indexes from [i*n, ((i+1)*n-1) + (p-1)]
        */
    
        pattern_value;
        hash_value = 1;
        
        int chars;                      /* number of characters that each row gets*/
        int i;                          /* loop index*/
        for(i = 0; i < p - 1; i++){
            chars = (i+1)*n + pattern_length - 2 - i*n;
            fread(elements, sizeof(char), chars, file);
            
            MPI_Send(elements, chars, MPI_CHAR, i, 1, MPI_COMM_WORLD);
        }

        /* gets its own elements*/
        chars = file_size - i*n;
        fread(elements, sizeof(char), chars, file);

        // pow(base, pattern_length-1)%p so it doesn't overflow
        for(i = 0; i < pattern_length - 1; i++){
            hash_value = (hash_value * base) % prime;
        }

        // calculates hash value for pattern
        for(i = 0; i < pattern_length; i++){
            pattern_value = (pattern_value * base + pattern[i]) % prime;
        }

        fclose(file);
    } else {
        MPI_Recv(elements, chars_for_p(id, n), MPI_CHAR, ROOT, 1, MPI_COMM_WORLD, &status);
    }
    // printf("%s \n%s \n", pattern, file_name);
    // /* Print size of file . intmax_t is a type that holds big numbers */
    // printf ("%8jd \n", (intmax_t)statbuf.st_size);
    int count = strlen(elements) - (pattern_length+1);
    char found[count];
    search(pattern, elements, pattern_length, pattern_value, hash_value);

    int prompt ; /* synchronizing variable */
    const int PROMPT_MSG = 2;
    const int RESPONSE_MSG = 3;
    if(0 == id){
        print_found();
        if(p > 1){
            int i;
            for(i = 1; i < p; i++){
                MPI_Send(&prompt, 1, MPI_INT, i, PROMPT_MSG, MPI_COMM_WORLD);

                MPI_Recv(found, count, MPI_CHAR, i, RESPONSE_MSG, MPI_COMM_WORLD, &status);
                print_found();
            }
        }
    } else {
        MPI_Recv(&prompt, 1, MPI_INT, 0, PROMPT_MSG, MPI_COMM_WORLD, &status);
        MPI_Send(found, count, MPI_CHAR, 0, RESPONSE_MSG, MPI_COMM_WORLD);
    }

    MPI_Finalize(); 
    return 0;
}

void search(char* pattern, char* text, int pattern_length, int pattern_value, int hash_value){
    int text_value = 0; /* Text hash value */

    int i, j;
    for(i = 0; i < pattern_length; i++){
        text_value = (text_value * base + text[i]) % prime;
    }

    for(i = 0; i <= strlen(text) - pattern_length; i++){
        if (pattern_value == text_value) {
            bool matches = true;
            for (j = 0; j < pattern_length; j++) {
                if (text[i + j] != pattern[j]){
                    matches = false;
                    break;
                }
            }
            if (matches) printf("%d \n", i);
        }
        text_value = (base * (text_value - text[i] * hash_value) + text[i + pattern_length]) % prime;
        
        /* if t can be negative*/
        if (text_value < 0) text_value += prime;
    }
}