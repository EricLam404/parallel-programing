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

#define ROOT 0
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
    int pattern_length;
    int pattern_value;
    int hash_value;

    MPI_Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &id );
    MPI_Comm_size (MPI_COMM_WORLD, &p);

    MPI_Bcast(&pattern_length, 1, MPI_INT, ROOT, MPI_COMM_WORLD); 
    MPI_Bcast(&pattern_value, 1, MPI_INT, ROOT, MPI_COMM_WORLD); 
    MPI_Bcast(&hash_value, 1, MPI_INT, ROOT, MPI_COMM_WORLD); 
    MPI_Bcast(&pattern_length, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&file_size, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(pattern, pattern_length, MPI_CHAR, ROOT, MPI_COMM_WORLD);
    MPI_Barrier (MPI_COMM_WORLD);
    
    /* Check usage */
    if (argc != 3) {
        if(ROOT == id){
            printf("Usage: %s <pattern> <file>\n", argv[0]);
        }
        MPI_Finalize(); 
        exit(1);
    }

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

    if (ROOT == 0) {
        for(int i = 0; i < file_size; i++){
            text[i] = fgetc(file);
        }
        pattern_length = strlen(pattern);
        pattern_value;
        hash_value = 1;
        int i;

        // pow(base, pattern_length-1)%p so it doesn't overflow
        for(i = 0; i < pattern_length - 1; i++){
            hash_value = (hash_value * base) % prime;
        }

        // calculates hash value for pattern
        for(i = 0; i < pattern_length; i++){
            pattern_value = (pattern_value * base + pattern[i]) % prime;
        }
    }
    
    // printf("%s \n%s \n", pattern, file_name);
    // /* Print size of file . intmax_t is a type that holds big numbers */
    // printf ("%8jd \n", (intmax_t)statbuf.st_size);
    search(pattern, text, pattern_length, pattern_value, hash_value);

    fclose(file);
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
