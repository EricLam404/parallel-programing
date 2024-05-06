/*******************************************************************************
  Title          : redact.c
  Author         : Eric Lam
  Created on     : May 1, 2024
  Description    : Redacting parrts of a text based on a pattern
  Purpose        : To practice how to use pthreads and matrixes
  Usage          : redact <number of threads> <pattern> <input file> <output file>
  Build with     : gcc -Wall -o redact redact.c -pthread
  Modifications  : May 6, 2024
  Added Error handling

*******************************************************************************/

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdint.h>
#include <stdbool.h>

#define SUCCESS             0
#define MALLOC_ERROR        1

char *text;
char *pattern;
int pattern_length;
pthread_mutex_t mutex;
pthread_barrier_t barrier;
int  **shared_region;
int  *shared_region_storage;

typedef struct _task_data
{
    int first;              /* index of first element in array for this thread */
    int last;               /* index of last element in array for this thread */
    int *array;             /* pointer to start of data array */
    int task_id;            /* program’s id of thread */
    pthread_t thread_id;    /* system’s id for thread */
    int result;             /* whatever the thread computes */
} task_data;

/**
 * Returns the redaction character based on the threads id
*/
char get_redaction_character(int rank) {
    char* redaction_array = "0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ_ ";
    return redaction_array[rank % 64];
}

/**
 * Finds all the indexes where the pattern matches the text, then redacts the indexes
 * based on the threads id
*/
void *redact_text(void  *thread_data);

/******************************************************************************/
void alloc_matrix(
        int     nrows,          /* number of rows in matrix                   */
        int     ncols,          /* number of columns in matrix                */
        size_t  element_size,   /* number of bytes per matrix element         */
        void  **matrix_storage, /* address of linear storage array for matrix */
        void ***matrix,         /* address of start of matrix                 */
        int    *errvalue)       /* return code for error, if any              */
{
    int   i;
    void *ptr_to_row_in_storage; /* pointer to a place in linear storage array
                                    where a row begins                        */
    void **matrix_row_start;     /* address of a 2D matrix row start pointer
                                    e.g., address of (*matrix)[row]           */
    size_t total_size;           /* amount of memory to allocate              */

    //printf("alloc_matrix called with r=%d,c=%d,e=%d\n",nrows, ncols, element_size);

    total_size = nrows * ncols * element_size;

    /* Step 1: Allocate an array of nrows * ncols * element_size bytes  */
    *matrix_storage = calloc(total_size, element_size);
    if ( NULL == *matrix_storage ) {
        /* malloc failed, so set error code and quit */
        *errvalue = MALLOC_ERROR;
        return;
    }

    /* Step 2: To create the 2D matrix, first allocate an array of nrows
       void* pointers */
    *matrix = malloc (nrows * sizeof(void*));
    if ( NULL == *matrix ) {
        /* malloc failed, so set error code and quit */
        *errvalue = MALLOC_ERROR;
        return;
    }


    /* Step 3: (The hard part) We need to put the addresses into the
       pointers of the 2D matrix that correspond to the starts of rows
       in the linear storage array. The offset of each row in linear storage
       is a multiple of (ncols * element_size) bytes.  So we initialize
       ptr_to_row_in_storage to the start of the linear storage array and
       add (ncols * element_size) for each new row start.
       The pointers in the array of pointers to rows are of type void*
       so an increment operation on one of them advances it to the next pointer.
       Therefore, we can initialize matrix_row_start to the start of the
       array of pointers, and auto-increment it to advance it.
    */

    /* Get address of start of array of pointers to linear storage,
       which is the address of first pointer, (*matrix)[0]   */
    matrix_row_start = (void*) &(*matrix[0]);

    /* Get address of start of linear storage array */
    ptr_to_row_in_storage = (void*) *matrix_storage;

    /* For each matrix pointer, *matrix[i], i = 0... nrows-1,
       set it to the start of the ith row in linear storage */
    for ( i = 0; i < nrows; i++ ) {
        /* matrix_row_start is the address of (*matrix)[i] and
           ptr_to_row_in_storage is the address of the start of the
           ith row in linear storage.
           Therefore, the following assignment changes the contents of
           (*matrix)[i]  to store the start of the ith row in linear storage
        */
        *matrix_row_start = (void*) ptr_to_row_in_storage;

        /* advance both pointers */
        matrix_row_start++;     /* next pointer in 2d array */
        ptr_to_row_in_storage +=  ncols * element_size; /* next row */
    }
    *errvalue = SUCCESS;


}

int main (int argc, char *argv[])
{
    int             retval;             /* return value of pthread_create   */
    int             t;                  /* thread loop variable             */
    struct stat     statbuf;            /* statbuf of the text file         */
    int             file_size;          /* size of the entire file in chars */
    FILE            *input_file;        /* input file location              */
    FILE            *output_file;       /* output file location             */
    int             num_threads;        /* total number of threads          */
    char            *input_file_name;   /* name of input file               */
    char            *output_file_name;  /* name of output file              */
    int             error;              /* error variable                   */

    task_data       *thread_data;  /* dynamically allocated array of thread data */
    pthread_attr_t  attr;

    char c;
    int temp;

    /* Make all threads joinable */
    pthread_attr_init(&attr);

    /* Check usage and errors */
    if (argc != 5) {
        printf("Usage: %s <number of threads> <pattern> <input file> <output file>\n", argv[0]);
        exit(1);
    }

    char *endptr;
    num_threads = (int)strtol(argv[1], &endptr, 10);
    if(*endptr != '\0' || num_threads <= 0){
        printf("%s is not a postive integer\n", argv[2]);
        exit(1);
    }

    pattern = argv[2];
    pattern_length = strlen(pattern);
    input_file_name = argv[3];
    output_file_name = argv[4];
    
    /* Call stat () to fill statbuf struct */ 
    if (lstat(input_file_name, &statbuf ) == -1 ) {
        printf("Could not stat file %s \n" , input_file_name);
        exit(1);
    }

    file_size = statbuf.st_size;
    text = calloc(file_size, sizeof(char));
    
    input_file = fopen(input_file_name, "r");
    temp = 0;
    while ((c = fgetc(input_file)) != EOF) {
        text[temp] = c;
        temp++;
    }

    output_file = fopen(output_file_name, "w");
    if (output_file == NULL) {
        perror("Error opening output file");
        exit(1);
    }

    /* Allocate the array of task_data structures on the heap.
        This is necessary because the array is not global.   */
    thread_data = calloc( num_threads, sizeof(task_data));

    if  ( thread_data == NULL  )
        exit(1);

    /* allocates the region where the threads share indexes*/
    alloc_matrix( num_threads, pattern_length-1, sizeof(int),
                 (void**) &shared_region_storage, (void***) &shared_region, &error);

    pthread_mutex_init(&mutex, NULL);
    /* Initialize a barrier with a count equal to the number of threads */
    pthread_barrier_init(&barrier, NULL, num_threads);

    /* Initialize task_data for each thread and then create the thread */
    int n = file_size - pattern_length + 1;
    for ( t = 0 ; t < num_threads; t++) {
        thread_data[t].first     = (t*n)/num_threads;
        thread_data[t].last      = ((t+1)*n)/num_threads - 1;
        thread_data[t].task_id   = t;
        retval = pthread_create(&(thread_data[t].thread_id), &attr,
                                  redact_text, (void *) &thread_data[t]);
        if ( retval ) {
            fprintf(stderr, "ERROR; return code from pthread_create() is %d\n", retval);
            exit(-1);
        }
    }

   /* Join all threads  */
    for ( t = 0 ; t < num_threads; t++) {
        pthread_join(thread_data[t].thread_id, (void**) NULL);
    }

    fputs(text, output_file);
    /* Free all memory allocated to program */
    free ( thread_data );
    fclose(input_file);
    fclose(output_file);
    return 0;   
}

// 1234567890 12345
// 12 34 56
// 123456  345678
//   [4] [4]
void *redact_text(void  *thread_data) {
    task_data* data = (task_data*) thread_data;
    char redaction_char = get_redaction_character(data->task_id);
    int i, j;
    int retval;

    /* dynamic array for matching indexes*/
    int size = 0;
    int *matching_indexes = calloc((data->last - data->first + 1), sizeof(int));

    // printf("Thread %d: (%d, %d) ", data->task_id, data->first, data->last);
    // for(i = data->first; i <= data->last; i++){
    //     printf("%c", text[i]);
    // }
    // printf("\n");

    /* loops through the threads portion of the text and finds the matching indexes*/
    for (i = data->first; i <= data->last; i++){
        for(j = 0; j < pattern_length; j++){
            if(text[i + j] != pattern[j]){
                break;
            }
        }

        if(j == pattern_length){
            matching_indexes[size] = i;
            size++;
        }
    }

    /* puts a barrier to wait for every thread to finish*/
    retval = pthread_barrier_wait(&barrier);
    if ( PTHREAD_BARRIER_SERIAL_THREAD != retval && 0 != retval )
        pthread_exit((void*) 0);

    /* loops through and redacts the text*/
    int index;
    for(i = 0; i < size; i++){
        index = matching_indexes[i];
        if(index < data->first+pattern_length || index > data->last-pattern_length){
            pthread_mutex_lock(&mutex);
            for(j = 0; j < pattern_length; j++){
                if(index + j > data->last){
                    if(shared_region[data->task_id][index+j-data->last-1] > data->task_id){
                        break;
                    } else {
                        shared_region[data->task_id][j] = data->task_id;
                    }
                }
                text[index + j] = redaction_char;
                if(index+j < data->first+pattern_length && data->task_id != 0){
                    shared_region[data->task_id-1][index+j - data->first] = data->task_id;
                }
            }
            pthread_mutex_unlock(&mutex);
        } else {
            for(j = 0; j < pattern_length; j++){
                text[index + j] = redaction_char;
            }
        }
    }

    pthread_exit((void*) 0);
}
