/******************************************************************************
Title          : steady_state.c
Author         : Eric Lam
Created on     : April  10, 2024
Description    : Uses a random walk to find the steady-state heat distribution
Purpose        : To demonstrate how to perform a random walk using a Dirichlet
                problem as input with parallelization.
Usage          : steady_state <inputfile> <x-cord> <y-cord>
                where <inputfile> is a plain text file containing:
                the dimensions of a 2d matrix, followed by 4 floats,
                representing the temperatures on the N, E, S, and W edges
                of the plate 
                <x-cord> <y-cord> are integers representing the cordinate of
                the point on the matrix to print the temperature of
Build with     : mpicc -g -Wall -o steady_state steady_state.c
Modifications  : April 14, 2024
                Fixed syntax issues where I was assigning values to a variable
                instead of comparing
******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <errno.h>
#include <string.h>
#include "mpi.h"

#define SUCCESS             0
#define MALLOC_ERROR        1
/*
Input file consists of
integer N  size of square array
four doubles: temperatures at N, E, S, and W in that order
*/

#define CONVERGENCE_THRESHOLD  0.05

#define NORTH   1
#define EAST    2
#define SOUTH   3
#define WEST    4

#define ROOT    0

typedef struct point2d_tag {
    int x;
    int y;
} point2d;

const point2d East  = {1, 0};
const point2d West  = {-1,0};
const point2d North = {0, 1};
const point2d South = {0,-1};

double uniform_random()
{
    return (double) (random()) / RAND_MAX;
}

point2d next_dir()
{
    double u;

    u = uniform_random();
    if ( u < 0.25 )
        return North;
    else if ( u < 0.50 )
        return East;
    else if ( u < 0.75 )
        return South;
    else
        return West;
}

int  on_boundary(point2d point, int width, int height)
{
    if ( 0 == point.x )
        return WEST;
    else if ( width -1 == point.x )
        return EAST;
    else if ( 0 == point.y )
        return NORTH;
    else if ( height - 1 == point.y )
        return SOUTH;
    else
        return 0;
}

point2d next_point(point2d oldpoint, point2d direction)
{
    point2d temp;

    temp.x = oldpoint.x + direction.x;
    temp.y = oldpoint.y + direction.y;
    return temp;
}

/** init_random()  initializes the state for the C random() function
 *  @param  int    state_size  [IN]  Size of state array for random to use
 *  @return char*  a pointer to the state array allocated for random()
 *  @post          After this call, an array of size state_size*sizeof(char) has
 *                 been allocated and initialized by C initstate(). It must be
 *                 freed by calling free()
 */
char*  init_random( int state_size, int id )
{
    char * state;
    const int seed = 1713050873;
    state  = (char*) malloc ( state_size * sizeof(char));
    if ( NULL != state )
        initstate((seed + id), state, state_size);
    return state;
}

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

int main ( int argc, char * argv[] )
{
    point2d current, next;      /* points for the random walk                   */
    int     i, j,  error;       /* looping variables and error                  */
    int     width, height;      /* width and height of the entire matrix        */
    double  oldvalue, diff;     /* old average and difference                   */
    double  maxdiff;            /* max difference in one loop                   */
    double  tolerance;          /* tolerance before reaching a stable state     */
    char    *randtable;         /* table to hold the inital values of the RNG   */
    int     tablesize = 1000;   /* RNG table size                               */
    FILE    *inputfile;         /* file with the input                          */
    int     location;           /* location of the current point                */
    double  boundary_temperature[4]; /* NESW boundary temperature               */
    double  **plate;            /* matrix of the points                         */
    double  *plate_storage;     /* storage of the matrix                        */
    int     id;                 /* process rank in MPI_COMM_WORLD               */
    int     p;                  /* number of processes in MPI_COMM_WORLD        */
    int     x,y;                /* cords of the point to find                   */
    int     local_height;       /* height of the local matrix                   */
    int     local_height_max;   /* max height in the local processor to check   */
    MPI_Status status;          /* result of MPI_Recv                           */

    MPI_Init (&argc, &argv);
    MPI_Comm_rank (MPI_COMM_WORLD, &id);
    MPI_Comm_size (MPI_COMM_WORLD, &p);

    // checks if there are any incorrect usages and handles the error
    if(ROOT == id){
        if ( argc < 4 ) {
            fprintf(stderr, "Usage: %s file xcord ycord\n", argv[0]);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        char *endptr;
        y = (int)strtol(argv[2], &endptr, 10);
        if(*endptr != '\0'){
            fprintf(stderr, "%s is not a integer\n", argv[2]);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        x = (int)strtol(argv[3], &endptr, 10);
        if(*endptr != '\0'){
            fprintf(stderr, "%s is not a integer\n", argv[3]);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        inputfile = fopen( argv[1],  "r");
        if (NULL == inputfile ) {
            fprintf(stderr, "Failed to open %s\n", argv[1] );
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        
        fscanf(inputfile, "%d ", &height);
        fscanf(inputfile, "%d ", &width);

        if(!(0 <= y && y < height)){
            fprintf(stderr, "%d is not a in between [0, %d)\n", y, height);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        if(!(0 <= x && x < width)){
            fprintf(stderr, "%d is not a in between [0, %d)\n", x, width);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        if(height < p){
            fprintf(stderr, "%d is less than the number of processors\n", height);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        for ( i = 0; i < 4; i++ )
            fscanf(inputfile, " %lf ", &(boundary_temperature[i]) ); 
        
        fclose(inputfile);
    }

    // broadcasts the needed variables
    MPI_Bcast(&width, 1, MPI_INT, ROOT, MPI_COMM_WORLD); 
    MPI_Bcast(&height, 1, MPI_INT, ROOT, MPI_COMM_WORLD); 
    MPI_Bcast(&x, 1, MPI_INT, ROOT, MPI_COMM_WORLD); 
    MPI_Bcast(&y, 1, MPI_INT, ROOT, MPI_COMM_WORLD); 
    MPI_Bcast(boundary_temperature, 4, MPI_DOUBLE, ROOT, MPI_COMM_WORLD); 

    // each processor gets their own local rows
    int num_rows = height / p;
    if(id != p - 1){
        local_height = num_rows;
    } else {
        local_height = height - (num_rows * (p-1));
    }


    // allocates memory for the matrix
    alloc_matrix( local_height, width, sizeof(double),
                 (void**) &plate_storage, (void***) &plate, &error);
    if ( error != 0 ){
        fprintf(stderr, "Alloc matrix error \n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    randtable = init_random(tablesize, id);
    /* Initialize temperatures at the four corners to be the average of the
       temperatures of the adjacent edges in the different processors
       and sets the NORTH and SOUTH boundaries to the input temp
    */
    if(id == 0){
        plate[0][0]              = (boundary_temperature[0] + boundary_temperature[3])/2;
        plate[0][width-1]        = (boundary_temperature[0] + boundary_temperature[1])/2;
        for ( j = 1; j < width - 1; j++ ) {
            plate[0][j]        = boundary_temperature[0];
        }
    }
    if(id == p - 1){
        plate[local_height-1][0]       = (boundary_temperature[3] + boundary_temperature[2])/2;
        plate[local_height-1][width-1] = (boundary_temperature[2] + boundary_temperature[1])/2;
        for ( j = 1; j < width - 1; j++ ) {
            plate[local_height-1][j] = boundary_temperature[2];
        }
    }

    // sets the EAST AND WEST boundary temps in every processor
    i = 0;
    local_height_max = local_height;
    if(id == 0) i = 1;
    else if(id == p - 1) local_height_max = local_height - 1;

    for (; i < local_height_max; i++ ) {
        plate[i][0]        = boundary_temperature[3];
        plate[i][width-1]  = boundary_temperature[1];
    }


    // does the random walk in each processors matrix until
    // a steady state is reached
    int count = 0;
    tolerance = CONVERGENCE_THRESHOLD;
    while ( 1 ) {
        maxdiff = 0;

        i = 0;
        local_height_max = local_height;
        if(id == 0) i = 1;
        else if(id == p - 1) local_height_max = local_height - 1;

        for (; i < local_height_max; i++ ) {
            for ( j = 1; j < width - 1 ; j++ ) {
                current.x = j;
                current.y = i + id*num_rows;
                while ( 0 == (location = on_boundary(current, width, height)) ) {
                    next = next_point(current, next_dir());
                    current = next;
                }
                oldvalue = plate[i][j];
                plate[i][j] = (oldvalue*count + boundary_temperature[location-1]) /
                              (count + 1);
                diff = fabs(plate[i][j] - oldvalue);
                if ( diff > maxdiff )
                    maxdiff = diff;
                /* maxdiff is largest difference in current iteration */
            }
        }
        if ( maxdiff < tolerance ){
            break;
        }
        else {
            count++;
        }
    }

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
    
    int prompt; /* synchronizing variable */
    const int PROMPT_MSG = 2;
    const int RESPONSE_MSG = 3;

    // finds the correct processor that has the cords of the point
    int p_send = y/num_rows;
    if(p_send > p-1) p_send = p-1;
    // prints if the processor is the ROOT else it sends a mpi 
    // signal to the processor that has the point
    if(ROOT == p_send && ROOT == id){
        printf("%.2f\n", plate[y][x]);
    } else{
        if(ROOT == id){
            double temp;
            MPI_Send(&prompt, 1, MPI_INT, p_send, PROMPT_MSG, MPI_COMM_WORLD);
            MPI_Recv(&temp, 1, MPI_DOUBLE, p_send, RESPONSE_MSG, MPI_COMM_WORLD, &status);
            printf("%.2f\n", temp);
        } else if(p_send == id) {
            /* waits for the processors turn to send data to the 0 processor */
            double temp = plate[y - (num_rows * id)][x];
            MPI_Recv(&prompt, 1, MPI_INT, ROOT, PROMPT_MSG, MPI_COMM_WORLD, &status);
            MPI_Send(&temp, 1, MPI_DOUBLE, ROOT, RESPONSE_MSG, MPI_COMM_WORLD);
        }
    }
    
    // printing out matrix for tesitng purposes
    // const int PRINT_MSG = 4;
    // const int FINISHED_MSG = 5;

    // if(0 == id){
    //     for ( i = 0; i < local_height; i++ ) {
    //         printf("pid:%d\t", id);
    //         for ( j = 0; j < width; j++ ) {
    //             printf("%.2f\t", plate[i][j]);
    //         }
    //         printf("\n");
    //     }

    //     for(i = 1; i < p; i++){
    //         MPI_Send(&prompt, 1, MPI_INT, i, PRINT_MSG, MPI_COMM_WORLD);
    //         MPI_Recv(&prompt, 1, MPI_INT, i, FINISHED_MSG, MPI_COMM_WORLD, &status);
    //     }
    // } else {
    //     /* waits for the processors turn to print matrix data to the 0 processor */
    //     MPI_Recv(&prompt, 1, MPI_INT, 0, PRINT_MSG, MPI_COMM_WORLD, &status);
    //     for ( i = 0; i < local_height; i++ ) {
    //         printf("pid:%d\t", id);
    //         for ( j = 0; j < width; j++ ) {
    //             printf("%.2f\t", plate[i][j]);
    //         }
    //         printf("\n");
    //     }
    //     MPI_Send(&prompt, 1, MPI_INT, 0, FINISHED_MSG, MPI_COMM_WORLD);

    // }

    free( plate);
    free( plate_storage );
    free(randtable);
    MPI_Finalize();
    return 0;
}
