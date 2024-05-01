/******************************************************************************
Title : search.c
Author : Eric Lam
Created on : March 24, 2024
Description : Estimates the natural log of a number based on the number of 
              segments
Purpose : Demonstrates parallel programing using MPI by splitting the fixe file
          into multiple text fragments and passing it through to each processor
          calling MPI_Send and MPI_RECV to gather the results back to the ROOT
          processor to print the results
Usage : finds the indexes of pattern matches in the text file
Build with : mpicc -Wall -g -o search search.c
Modifications: March 27, 2024
Fixed malloc error with larger text files
******************************************************************************/

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
#include <ctype.h>
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
int* search(
    char* pattern, 
    char* text, 
    int pattern_length, 
    int pattern_value, 
    int hash_value, 
    int *size
);

/**
 * prints the indexes in the array where a match was found
 * n is the offset so that it prints the original index
 * found in file
*/
void print(int *arr, int size, int id, int n);

/**
 * checks if the string contains all non control characters
*/
bool isValid(char *string, int size);

int main(int argc, char *argv[]) {
    struct stat statbuf;    /* statbuf of the text file     */
    int id;                 /* rank of executing process    */
    int p;                  /* number of processes          */
    int pattern_length;     /* length of pattern            */
    int pattern_value;      /* hash value of pattern        */
    int hash_value;         /* general hash value           */
    int n;                  /* number of chars for each p   */

    int file_size;         
    FILE *file;
    MPI_Status status;      /* result of MPI_Recv */

    MPI_Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &id );
    MPI_Comm_size (MPI_COMM_WORLD, &p);
    
    const int ROOT = p - 1; /* root processor is the last processor*/
    int i;                  /* loop index*/
    if(ROOT == id){
        if((argc < 3)){
            printf("Usage: %s <pattern> <file>\n", argv[0]);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        /* Call stat () to fill statbuf struct */ 
        if (lstat(argv[2], &statbuf ) == -1 ) {
            printf("Could not stat file %s \n" , argv[2]);
            MPI_Abort(MPI_COMM_WORLD, 1);
        } 
        file_size = statbuf.st_size;
        pattern_length = strlen(argv[1]);
        if(!isValid(argv[1], pattern_length)){
            printf("The pattern contains non-valid characters \n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        if(file_size < pattern_length){
            printf("File size is less than pattern length \n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        /* Rounding n down to the nearest int*/
        n = (((double)file_size-(pattern_length-1)) / p); 
    }
    
    MPI_Bcast(&pattern_length, 1, MPI_INT, ROOT, MPI_COMM_WORLD); 
    char *pattern = (char *)malloc((pattern_length + 1) * sizeof(char));
    if (!pattern) {
        printf("Memory allocation failed.\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    if(ROOT == id){
        strcpy(pattern, argv[1]);
        hash_value = 1;
        pattern_value = 0;

        // pow(base, pattern_length-1)%p so it doesn't overflow
        for(i = 0; i < pattern_length - 1; i++){
            hash_value = (hash_value * base) % prime;
        }

        // calculates hash value for pattern
        for(i = 0; i < pattern_length; i++){
            pattern_value = (pattern_value * base + pattern[i]) % prime;
        }
    }
    
    MPI_Bcast(pattern, pattern_length + 1, MPI_CHAR, ROOT, MPI_COMM_WORLD); 

    MPI_Bcast(&hash_value, 1, MPI_INT, ROOT, MPI_COMM_WORLD); 
    MPI_Bcast(&pattern_value, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&n, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&file_size, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
    /**
     * if the number of indexes to check is smaller than the number of processors
     * only the root processor performs the search
    */
    if(file_size - (pattern_length-1) < p){
        if(ROOT == id){
            file = fopen(argv[2], "r");

            char *elements = (char *)calloc(file_size, sizeof(char));
            if (!elements) {
                printf("Memory allocation failed.\n");
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
            fread(elements, sizeof(char), file_size, file);

            int size = 0;
            int *match = search(
                            pattern, 
                            elements, 
                            pattern_length, 
                            pattern_value, 
                            hash_value, 
                            &size);
            print(match, size, 0, n);
        }   

        free(pattern);
        MPI_Finalize();
        return 0;
    }

    char *elements = (char *)calloc((file_size-((p-1)*n)), sizeof(char));
    if (!elements) {
        printf("Memory allocation failed.\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    int chars;                      /* number of characters that each row gets*/
    /* elements to send*/
     #ifdef DEBUG_ON
    /* To debug, compile this program with the -DDEBUG_ON option,
    which defines the symbol DEBUG_ON, and run the program as usual
    with mpirun.
    When the output appears on the terminal, listing the pids of the
    processes and which hosts they are on, choose the lowest
    pid P on the machien you are connected to.
    Open a new terminal window and in that window issue the command
    gdb --pid P
    (or gdb -p P on some systems)
    and after gdb starts, go up the stack to main by entering the
    command
    up 3
    (main will be three stacks frames above your current frame,
    which should be nanosleep.)
    Then enter the command
    set var i = 1
    to break the while loop. You can now run ordinary gdb commands to
    debug this process. This should be process 0.
    Repeat these steps for each other process that you created in
    the mpirun command.
    */
    #include <unistd.h>
    int debug = 0;
    char hostname[256];
    gethostname(hostname, sizeof(hostname));
    printf("PID %d on %s ready for attach\n", getpid(), hostname);
    fflush(stdout);
    while (0 == debug)
    sleep(5);
    #endif
    chars = n+pattern_length-1;
    if(ROOT == id){
        file = fopen(argv[2], "r");

        /**
         * send the indexes from [0, n+pattern_length-1]
        */
        fread(elements, sizeof(char), chars, file);
        MPI_Send(elements, chars, MPI_CHAR, 0, 1, MPI_COMM_WORLD);


        /**
         * each processor gets [n*i, (n+pattern_length-1) * i + 1]
         * each prosssor only moves foward n characters, so that
         * there are enough characters for each processor to check
         * n times
        */
        
        char *temp = (char *)calloc(chars, sizeof(char));
        for(i = 1; i < p - 1; i++){
            strcpy(temp, elements+n);
            fread(elements, sizeof(char), pattern_length-1, file);
            strcat(temp, elements);
            MPI_Send(temp, chars, MPI_CHAR, i, 1, MPI_COMM_WORLD);
        }

        /* gets its own elements*/
        chars = file_size - i*n;
        strcpy(temp, elements+n);
        free(elements);
        elements = (char *)calloc(chars, sizeof(char));
        if (!elements) {
            printf("Memory allocation failed.\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        fread(elements, sizeof(char), chars, file);
        strcat(temp, elements);
        strcpy(elements, temp);
        free(temp);
        fclose(file);
    } else {
        MPI_Recv(elements, chars, MPI_CHAR, ROOT, 1, MPI_COMM_WORLD, &status);
    }

    int size = 0;
    int *match = search(
                    pattern, 
                    elements, 
                    pattern_length, 
                    pattern_value, 
                    hash_value, 
                    &size
                );
    int prompt; /* synchronizing variable */
    const int PROMPT_MSG = 2;
    const int SIZE_MSG = 3;
    const int RESPONSE_MSG = 4;

    if(0 == id){
        /* prints the outputs in order */
        print(match, size, id, n);
        for(i = 1; i < p; i++){
            MPI_Send(&prompt, 1, MPI_INT, i, PROMPT_MSG, MPI_COMM_WORLD);
            MPI_Recv(&size, 1, MPI_INT, i, SIZE_MSG, MPI_COMM_WORLD, &status);

            MPI_Recv(match, size, MPI_INT, i, RESPONSE_MSG, MPI_COMM_WORLD, &status);
            print(match, size, i, n);
        }
    } else {
        /* waits for the processors turn to send data to the 0 processor */
        MPI_Recv(&prompt, 1, MPI_INT, 0, PROMPT_MSG, MPI_COMM_WORLD, &status);
        MPI_Send(&size, 1, MPI_INT, 0, SIZE_MSG, MPI_COMM_WORLD);
        MPI_Send(match, size, MPI_INT, 0, RESPONSE_MSG, MPI_COMM_WORLD);
    }

    free(match);
    free(pattern);
    MPI_Finalize();
    return 0;
}

int* search(
    char* pattern, 
    char* text, 
    int pattern_length, 
    int pattern_value, 
    int hash_value, 
    int *size
){
    int text_value = 0;     /* Text hash value */
    int i, j;
    /* alloc enough memory for if every index is a match*/
    int* array = (int*)calloc(strlen(text) + (pattern_length-1), sizeof(int));
    if (!array) {
        printf("Memory allocation failed.\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    /* calculates the text index */
    for(i = 0; i < pattern_length; i++){
        text_value = (text_value * base + text[i]) % prime;
    }

    /**
     * finds potention string matches using a rolling hash
     * if the hash of the text is the same as the pattern 
     * then you perform a check to see if the pattern 
     * matches the string segment from 
     * point i to i + pattern_length
     * if it is then it adds it to the array and updates
     * size to contain the true number of elements
    */
    for(i = 0; i <= strlen(text) - pattern_length; i++){
        if (pattern_value == text_value) {
            bool matches = true;
            for (j = 0; j < pattern_length; j++) {
                if (text[i + j] != pattern[j]){
                    matches = false;
                    break;
                }
            }
            if (matches){
                array[*size] = i;
                (*size)++;
            }
        }
        text_value = (base * (text_value - text[i] * hash_value) 
                        + text[i + pattern_length]) % prime;
        
        /* if t can be negative*/
        if (text_value < 0) text_value += prime;
    }

    return array;
}

void print(int *arr, int size, int id, int n){
    int offset = id*n;
    for(int i = 0; i < size; i++){
        printf("%d \n", (offset+arr[i]));
    }
}

bool isValid(char *string, int size){
    for(int i = 0; i < size; i++){
        if (isprint(string[i]) == 0){
            return false;
        }
    }
    return true;
}