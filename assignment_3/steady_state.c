/******************************************************************************
  Title          : steady_state.c
  Author         : Eric Lam, modifered from Professor Stewart Weiss randomwalk
  Created on     : April  13, 2024
  Description    : Uses a random walk to find the steady-state heat distribution
  Purpose        : To demonstrate how to performa  random walk using a Dirichlet
                   problem as input.
  Usage          : steady_state <inputfile> <x-cord> <y-cord>
                   where <inputfile> is a plain text file containing:
                   the dimensions of matrix, followed by 4 floats,
                   representing the temperatures on the N, E, S, and W edges
                   of the plate
  Build with     : gcc -g -Wall -o steady_state steady_state.c
  Modifications  :

  Note:
  This will run for a very long time. For a 100 by 100 grid on a relatively
  fast processor, be prepared to wait 30 minutes or more. This is why it is
  better to use a parallel program to solve it.

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <errno.h>
#include <string.h>

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
char*  init_random( int state_size )
{
    char * state;
    state  = (char*) malloc ( state_size * sizeof(char));
    if ( NULL != state )
        initstate(time(NULL), state, state_size);
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
    point2d current, next;
    int     count = 0;
    int     i, j,  error;
    int     width, height;
    double  oldvalue, diff;
    double  maxdiff;
    double  tolerance;
    char    *randtable;
    int     tablesize = 1000;
    FILE    *inputfile;
    int     location;
    double  boundary_temperature[4];
    double  **plate;
    double  *plate_storage;
    int     x, y;

    if ( argc < 4 ) {
        fprintf(stderr, "Usage: %s file xcord ycord\n", argv[0]);
        exit(1);
    }

    char *endptr;
    x = (int)strtol(argv[2], &endptr, 10);
    if(*endptr != '\0'){
        fprintf(stderr, "%s is not a integer\n", argv[2]);
        exit(1);
    }

    y = (int)strtol(argv[3], &endptr, 10);
    if(*endptr != '\0'){
        fprintf(stderr, "%s is not a integer\n", argv[3]);
        exit(1); 
    }

    inputfile = fopen( argv[1],  "r");
    if (NULL == inputfile ) {
        fprintf(stderr, "Failed to open %s\n", argv[1] );
        exit(1);
    }
    
    fscanf(inputfile, "%d ", &width);
    fscanf(inputfile, "%d ", &height);

    if(!(0 <= x && x < width)){
        fprintf(stderr, "%d is not a in between [0, %d)\n", x, width);
        exit(1); 
    }
    if(!(0 <= y && y < height)){
        fprintf(stderr, "%d is not a in between [0, %d)\n", y, height);
        exit(1); 
    }

    for ( i = 0; i < 4; i++ )
        fscanf(inputfile, " %lf ", &(boundary_temperature[i]) ); 
    alloc_matrix( width, height, sizeof(double),
                 (void**) &plate_storage, (void***) &plate, &error);
    if ( error != 0 )
        exit(1);

    randtable = init_random(tablesize);
    /* Initialize temperatures at the four corners to be the average of the
       temperatures of the adjacent edges
    */
    plate[0][0]              = (boundary_temperature[0] + boundary_temperature[3])/2;
    plate[0][width-1]        = (boundary_temperature[0] + boundary_temperature[1])/2;
    plate[height-1][0]       = (boundary_temperature[3] + boundary_temperature[2])/2;
    plate[height-1][width-1] = (boundary_temperature[2] + boundary_temperature[1])/2;

    /* Initialize the termperatures along the edges of the plate */
    for ( j = 1; j < width -1; j++ ) {
        plate[0][j]        = boundary_temperature[0];
        plate[height-1][j] = boundary_temperature[2];
    }
    for ( i = 1; i < height -1; i++ ) {
        plate[i][0]        = boundary_temperature[3];
        plate[i][width-1]  = boundary_temperature[1];
    }

    count = 0;
    tolerance = CONVERGENCE_THRESHOLD;
    while ( 1 ) {
        maxdiff = 0;
        for ( i = 1; i < height -1 ; i++ ) {
            for ( j = 1; j < width -1 ; j++ ) {
                current.x = j;
                current.y = i;
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
        if ( maxdiff < tolerance )
            break;
        else {
            count++;
        }

        printf("\033[2J\033[H");
        printf("diff; %f < %f in %d rounds \n", maxdiff, tolerance, count);
    }
    printf("%.2f\n", plate[x][y]);

    free( plate);
    free( plate_storage );
    free(randtable);
    fclose(inputfile);
    return 0;
}

int size_of_block( int id, int ntotal_elements, int p )
{
    return ( ( ( id + 1) * ntotal_elements ) / p ) -
           ( ( id *      ntotal_elements ) / p );
}
