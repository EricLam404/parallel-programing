/******************************************************************************
  Title          : matrix_vector03.c
  Author         : Stewart Weiss, based on Michael Quinn's mv3.c
  Created on     : March  8, 2014
  Description    : Compute product of matrix A times vector X to standard output
  Purpose        : To show how to use a 2D block decomposition of a matrix
  Usage          : matrix_vector03  matrix-file  vector-file
  Build with     : mpicc -o matrix_vector03 -g -Wall -I../include -L../lib \
                         matrix_vector03.c -lutils
  Modifications  :
 
******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <string.h>

/* Change these two definitions when the matrix and vector
   element types change */

typedef int Element_type;
#define MPI_TYPE MPI_INT
#define REORDER  1
#define SUCCESS             0
#define  FLOAT  0
#define  DOUBLE 1
#define  INT    2
#define  CHAR   3
#define  LONG   4
#define  LLONG  5

//#define DEBUG_ON
void read_and_distribute_2dblock_matrix (
          char        *filename,       /* [IN]  name of file to read          */
          void      ***matrix,         /* [OUT] matrix to fill with data      */
          void       **matrix_storage, /* [OUT] linear storage for the matrix */
          MPI_Datatype dtype,          /* [IN]  matrix element type           */
          int         *nrows,          /* [OUT] number of rows in matrix      */
          int         *ncols,          /* [OUT] number of columns in matrix   */
          int         *errval,         /* [OUT] sucess/error code on return   */
          MPI_Comm     cart_comm);      /* [IN]  communicator handle           */
int size_of_block( int id, int ntotal_elements, int p );
int main (int argc, char *argv[]) 
{
    Element_type **A;         /* input matrix to be multiplied                */
    Element_type  *A_storage; /* backing storage for matrix                   */
    Element_type  *X_local;   /* subblock of input vector X in process        */
    Element_type  *B_sums;    /* subblock of result vector B = A*X            */
    Element_type  *B_partial; /* subblock of partial result vector            */
    Element_type *input_vec;  /* Dynamic storage for input vector in proc 0,0 */
    int       nrows;          /* total number of rows in input matrix         */
    int       ncols;          /* total number of columns in input matrix      */
    int       id;             /* process rank in MPI_COMM_WORLD               */
    int       p;              /* number of processes in MPI_COMM_WORLD        */
    int       i, j;           /* loop indices                                 */
    int       error;          /* error exit value of calls                    */
    int       nlocal_rows;    /* number of rows belonging to this process     */
    int       nlocal_cols;    /* number of cols belonging to this process     */
    int       grid_id;        /* rank of process within Cartesian grid        */
    int       row_id;         /* rank of process within its row subgroup      */
    int       grid_size[2];   /* grid dimensions                              */
    int       grid_coords[2]; /* process's coordinates in Cartesian grid      */
    MPI_Comm  grid_comm;      /* Cartesian communicator for grid of procs     */
    MPI_Comm  row_comm;       /* communicator for all processes in a grid row */
    MPI_Comm  col_comm;       /* communicator for all processes in a grid col */
    int      *send_count;     /* array of counts for MPI_Scatterv function    */
    int      *send_offset;    /* array of offsets for MPI_Scatterv function   */
    int       periodic[2];    /* array of wraparound flags when creating grid */
    int       nrows_X;        /* number of rows of X vector local to process  */
    FILE     *file;           /* for open input files                         */

    double    max_seconds;    /* for timing program              */
    double    seconds;        /* Elapsed time for matrix-vector multiply */

    MPI_Init (&argc, &argv);
    MPI_Comm_rank (MPI_COMM_WORLD, &id);
    MPI_Comm_size (MPI_COMM_WORLD, &p);

    /* Set up the various communicators */
    /* Let MPI_Dims_create choose the dimensions by setting grid_size[i]=0 */
    grid_size[0] = 0;
    grid_size[1] = 0;
    MPI_Dims_create (p, 2, grid_size);
    /* Now grid_size[0] will be number of rows in grid and grid_size[1] will
       be the number of columns in grid */

    /* The grid will not be periodic, so turn off wraparound flags */
    periodic[0]  = 0;
    periodic[1]  = 0;

    /* Create the Cartesian commincator from all processes in program. */
    MPI_Cart_create (MPI_COMM_WORLD, 2, grid_size, 
                     periodic, REORDER, &grid_comm);

    /* Get the rank within the Cartesian communicator of this process */
    MPI_Comm_rank (grid_comm, &grid_id);

    /* Get the coordinates within the grid of this process */
    MPI_Cart_coords (grid_comm, grid_id, 2, grid_coords);

    /* Each process now participates in a split of the Cartesian communicator
       into rows, preserving their relative ranking within each row. 
       When the third argument in the call to MPI_Comm_split is a unique
       value among all processes having the same second argument, that value
       is used as the process's rank in the newly created subgroup. Therefore,
       for all processes that supply the same grid row as argument 2, if they
       provide their column coordinate as argument 3, that column coordinate
       becomes their rank with the row subgroup. This implies that the process
       within each row in the first column will be process 0 within that new
       row group.
    */
    MPI_Comm_split (grid_comm, grid_coords[0], grid_coords[1], &row_comm);

    /* They do the same thing by columns now, preserving their relative
       rank within the column. The processes in the first row become the 
       processes with rank 0 in the column groups to which they belong. */
    MPI_Comm_split (grid_comm, grid_coords[1], grid_coords[0], &col_comm);

    /* All processes collectively acquire their sub-blocks of the matrix
       from the file. One process does the actual reading, and distributes
       to all other processes. This returns the number of rows and columns
       in the matrix through nrows and ncols. 
    */
    read_and_distribute_2dblock_matrix (argv[1], (void *) &A,
      (void *) &A_storage, MPI_TYPE, &nrows, &ncols, &error, grid_comm);
    if ( 0 != error )
        MPI_Abort (MPI_COMM_WORLD, MPI_ERR_IO);
   
    /* Each process obtains the number of rows and the number of columns
       it owns. The row count uses the formula floor(i*nrows/p)  where
           i = grid coordinate[0], nrows= total rows, p=grid dimension[0]
       and column count uses floor(j*ncols/p) where
           j = grid coordinate[1], ncols= total cols, p=grid dimension[1]
    */
    nlocal_rows = size_of_block(grid_coords[0], nrows, grid_size[0]);
    nlocal_cols = size_of_block(grid_coords[1], ncols, grid_size[1]);

    /* Process 0 in row 0 reads the vector file and stores the entire vector
       locally. */
    if ( 0 == grid_coords[0] ) {
        /* Get rank and make sure only process 0 does this. */
        MPI_Comm_rank(row_comm, &row_id);
        if ( 0 == row_id ) { /* we could have also just checked column coord. */
            /* Attempt to open the file and abort the program if unsuccessful */
            file = fopen (argv[2], "r");
            if ( NULL == file ) {
                MPI_Abort (MPI_COMM_WORLD, MPI_ERR_FILE );
            }
            else { 
                /* File was opened successfully. Read vector length into temp.*/
                int temp;
                fread (&temp, sizeof(int), 1, file);
                /* Make sure vector length matches matrix width. */
                if ( temp != ncols ) {
                    fprintf(stderr, "Vector size and matrix size unmatched.\n");
                    MPI_Abort(MPI_COMM_WORLD, MPI_ERR_IO);
                }

                /* Allocate storage for entire vector in this process */
                input_vec = (Element_type*) malloc(ncols*sizeof(Element_type)); 
                if ( NULL == input_vec ) {
                    fprintf(stderr, "Could not allocate more storage.\n");
                    MPI_Abort(MPI_COMM_WORLD, MPI_ERR_NO_MEM);
                }
                /* Read entire vector into input_vec */  
                fread (input_vec, sizeof(Element_type), ncols, file);
            }  
        }
    }
    /* Get the number of elements of the vector belonging to current process 
       and allocate a vector to store those elements. */
    nrows_X = size_of_block (grid_coords[1], ncols, grid_size[1] );
    X_local = (Element_type *) malloc (nrows_X * sizeof(Element_type));
    if (NULL == X_local ) {
        fprintf(stderr, "Could not allocate more storage.\n");
        MPI_Abort(MPI_COMM_WORLD, MPI_ERR_NO_MEM);
    }

   /* The vector is scattered by process with rank 0 in row 0 to all other 
      processes in row 0. First we need to construct the count and offset arrays
      for the MPI_Scatterv function, so we allocate the arrays and call
      our utility function to fill them.*/
    if ( 0 == grid_coords[0] ) {
        /* Get id in this group -- do not assume it is 0 */
        if (grid_size[1] > 1) {
            send_count  = (int *) malloc (grid_size[1] * sizeof(int));
            send_offset = (int *) malloc (grid_size[1] * sizeof(int));
            if ( NULL == send_count || NULL == send_offset ) {
                fprintf(stderr, "Could not allocate more storage.\n");
                MPI_Abort(MPI_COMM_WORLD, MPI_ERR_NO_MEM);
            }
            init_communication_arrays(
                         grid_size[1], ncols, send_count, send_offset);
            MPI_Scatterv (input_vec, send_count, send_offset, MPI_TYPE, X_local,
                          nrows_X, MPI_TYPE, 0, row_comm);
        } 
        else {
            /* If the grid has a single column, then there are no other
               processes and process 0 just needs to copy the vector into
               the array X. */
            for (i = 0; i < ncols; i++) X_local[i] = input_vec[i];
        }
    }

    /* Row 0 processess broadcast their subvectors to all processes in same 
       column */
    MPI_Bcast (X_local,nrows_X, MPI_TYPE, 0, col_comm);


    MPI_Barrier (MPI_COMM_WORLD);
    seconds = - MPI_Wtime();


   /* Every process allocates storage for its own block of the 
      partial sums of the inner product of AX.  It needs two such blocks;
      one for its local sums and a second to use in a global reduction
      to store the results. */
    B_partial  = (Element_type *) malloc (nlocal_rows * sizeof(Element_type));
    B_sums   = (Element_type *) malloc (nlocal_rows * sizeof(Element_type));
    if (NULL == B_sums || NULL == B_partial ) {
        fprintf(stderr, "Could not allocate more storage.\n");
        MPI_Abort(MPI_COMM_WORLD, MPI_ERR_NO_MEM);
    }

    /* Each process does its own local inner product of its sub-blocks of A
       and the partial vector X, storing the result in the partial vector B. */
    for (i = 0; i < nlocal_rows; i++) {
         B_partial[i] = 0.0;
         for (j = 0; j < nlocal_cols; j++) {
              B_partial[i] += A[i][j] * X_local[j];
         }
    }

    /* All process in each row perform a simultaneous reduction of nlocal_rows
       many partial sums, storing the global sums into the B_sums array in the 
       process with rank 0 within the row communicator, which will be the 
       processes in the first column. */
    MPI_Reduce(B_partial, B_sums, nlocal_rows, MPI_TYPE, MPI_SUM, 0, row_comm);

    /* All processes in the first column print their portions of the result 
       vector. These are the ones that just received the global sums. */
    if (grid_coords[1] == 0) {
        print_block_vector (B_sums, MPI_TYPE, nrows, col_comm);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    seconds += MPI_Wtime();
    max_seconds = seconds;
   // MPI_Reduce (&seconds, &max_seconds, 1, MPI_DOUBLE, MPI_MAX, 0,
    //  MPI_COMM_WORLD);

    if (0 == id) {
        printf ("MV5) N = %d, Processes = %d, Time = %12.6f sec,",
         ncols, p, max_seconds);
        printf ("Mflop = %6.2f\n", 2*ncols*ncols/(1000000.0*max_seconds));
    }
    MPI_Finalize();
    return 0;
}

void read_and_distribute_2dblock_matrix (
          char        *filename,       /* [IN]  name of file to read          */
          void      ***matrix,         /* [OUT] matrix to fill with data      */
          void       **matrix_storage, /* [OUT] linear storage for the matrix */
          MPI_Datatype dtype,          /* [IN]  matrix element type           */
          int         *nrows,          /* [OUT] number of rows in matrix      */
          int         *ncols,          /* [OUT] number of columns in matrix   */
          int         *errval,         /* [OUT] sucess/error code on return   */
          MPI_Comm     cart_comm)      /* [IN]  communicator handle           */
{
    int    i,j,k;           /* various loop index variables                  */
    int    grid_id;         /* process rank in the cartesian grid            */
    int    p;               /* number of processes in the cartesian grid     */
    size_t element_size;    /* number of bytes in matrix element type        */
    int    mpi_initialized; /* flag to check if MPI_Init was called already  */
    FILE   *file;           /* input file stream pointer                     */
    int    nlocal_rows;     /* number of rows that calling process "owns"    */
    int    nlocal_cols;     /* number of columns that calling process "owns" */
    MPI_Status   status;    /* result of MPI_Recv call                       */
    int    dest_id;         /* rank of receiving process in cartesian grid   */
    int    grid_coord[2];   /* process coordinates in the grid               */
    int    grid_periodic[2];/* flags indicating if grid wraps around         */
    int    grid_size[2];    /* dimensions of grid                            */
    void*  buffer;          /* address of temp location to store rows        */
    int    block_coord[2];  /* coordinates in grid of current block          */
    void*  source_address;  /* address of block to be sent                   */
    void*  dest_address;    /* location where block is to be received        */
    
    /* Make sure we are being called by a program that init-ed MPI */
    MPI_Initialized(&mpi_initialized);
    if ( !mpi_initialized ) {
       *errval = -1;
       return;
    }
  
    /* Get process rank in grid and the number of processes in group */
    MPI_Comm_rank (cart_comm, &grid_id);  
    MPI_Comm_size (cart_comm, &p);            
   
    /* Get the number of bytes in a matrix element */
    element_size = get_size (dtype);
    if ( element_size <= 0 ) {
       *errval = -1;
       return;
    }
       
    /* Process 0 opens the file and read the number of rows and columns */
    if ( 0 == grid_id ) { 
        /* Process 0 opens the binary file containing the matrix and 
           reads the first two numbers, which are the number of rows and 
           columns respectively. */
        file = fopen (filename, "r");
        if ( NULL == file ) {
            *nrows = 0;
            *ncols = 0;
        }
        else { /* successful open */
            fread (nrows, sizeof(int), 1, file);
            fread (ncols, sizeof(int), 1, file);
        }      
    }

    /* Process 0 broadcasts the numbers of rows to all other processes. */
    MPI_Bcast (nrows, 1, MPI_INT, 0, cart_comm);

    /* All processes check value of *nrows; if 0 it indicates failed open */
    if ( 0 == *nrows  ) {
       *errval = -1;
       return;
    }

    /* Process 0 broadcasts the numbers of columns to all other processes. 
       No need to check whether *ncols is zero. */
    MPI_Bcast (ncols, 1, MPI_INT, 0, cart_comm);

    /* All processes obtain the grid's topology so they can determine
       their block sizes. */
    MPI_Cart_get (cart_comm, 2, grid_size, grid_periodic, grid_coord);

    /* Each process sets nlocal_rows = the number of rows the process owns and
       local_cols to the number of columns it owns.
       The number of rows depends on the process's row coordinate, *nrows, 
       and the number of grid rows in total. This implements the formula
             blocks = floor((i+1)*n/p) - floor(i*n/p)
       where i is grid coordinate in given dimension, n is either total
       number of rows, or total number of columns, and p is the number of
       processes in grid in the given dimension.
    */
    nlocal_rows = size_of_block( grid_coord[0], *nrows, grid_size[0] );
    nlocal_cols = size_of_block( grid_coord[1], *ncols, grid_size[1] );

    /* Each process creates its linear storage and 2D matrix for accessing 
       the elements of its assigned rows. It needs storage for a 2D matrix
       of nlocal_rows by nlocal_cols elements. */
    alloc_matrix( nlocal_rows, nlocal_cols, element_size,
                  matrix_storage, 
                  matrix,
                  errval);
    if ( SUCCESS != *errval ) {
         MPI_Abort (cart_comm, *errval);
    }

   /* Grid process 0 reads in the matrix one row at a time
      and distributes each row among the MPI processes. The first step is
      to allocate storage for one row of the matrix. */
    if ( 0 == grid_id ) {
        buffer = malloc (*nrows * element_size);
        if ( buffer == NULL ) { 
            MPI_Abort (cart_comm, *errval);
        }
    }

    /* This is the read and distribute loop. Process 0 will read a row
       and send it to the processes that are supposed to have it. It needs to
       break it into blocks, with successive blocks going to the processes
       in the same grid row but successive grid columns. */

    for (i = 0; i < grid_size[0]; i++) { /* for each grid row */
        /* Set block_coord[0] to the current grid row index */
        block_coord[0] = i;

        /* For every matrix row that is part of this grid row */
        for (j = 0; j < size_of_block(i, *nrows, grid_size[0] ); j++) {

            /* Process 0 reads  a row of the matrix */
            if ( 0 == grid_id ) {
                fread (buffer, element_size, *ncols, file);
            }
            /* Every process executes this loop. For each grid column within
               the current grid row ... */
            for (k = 0; k < grid_size[1]; k++) {
                block_coord[1] = k;

                /* Determine the grid id of the process in grid position
                   [i,k].  This is the destination process. Its id is returned 
                   in dest_id. */
                MPI_Cart_rank (cart_comm, block_coord, &dest_id);

                /* Process 0 needs to determine the start of the block to
                   be sent to process dest_id. This is the start in a row
                   with *ncols elements assigned to kth process out of
                   grid_size[1] many processes in that row. */
                if ( 0 == grid_id ) {
                    source_address = buffer + 
                                 ( (k*(*ncols))/grid_size[1] ) * element_size;
                    /* The process has to make sure it does not try to send to
                       itself. If so it does a memory copy instead. */
                    if (0 == dest_id ) {
                        /* It is sending to itself */
                        dest_address = (*matrix)[j];
                        memcpy (dest_address, source_address, 
                               nlocal_cols * element_size);                  
                    } 
                    else {
                        /* It is sending to another process */
                        int blocksize = size_of_block(k,*ncols, grid_size[1]);
                        MPI_Send (source_address,blocksize, dtype,
                                  dest_id, 0, cart_comm);
                    }
                }
                else if (grid_id == dest_id) {
                         MPI_Recv ((*matrix)[j], nlocal_cols, dtype, 0,
                          0, cart_comm, &status);
                }
            } /* end for k */
        } /* end for j */
    } /* for i */

    if (grid_id == 0) 
        free (buffer);
}

int size_of_block( int id, int ntotal_elements, int p )
{
    return ( ( ( id + 1) * ntotal_elements ) / p ) -
           ( ( id *      ntotal_elements ) / p );
}

/******************************************************************************/
/** get_type() returns a constant representing MPI_Datatype argument
 *  @param    MPI_Datatype  t
 *  @return   integer constant
 */
int get_type (MPI_Datatype t)
{
    if ( ( t == MPI_BYTE) ||
         ( t == MPI_CHAR) ||
         ( t == MPI_SIGNED_CHAR) ||
         ( t == MPI_UNSIGNED_CHAR) )
        return CHAR;

    if ( t == MPI_DOUBLE )
        return DOUBLE;
    if (t == MPI_FLOAT)
        return FLOAT;

    if ( ( t == MPI_INT) ||
         ( t == MPI_UNSIGNED) )
        return INT;

    if ( ( t == MPI_LONG) ||
         ( t == MPI_UNSIGNED_LONG) )
        return LONG;

    if ( ( t == MPI_LONG_LONG) ||
         ( t == MPI_LONG_LONG_INT) ||
         ( t == MPI_UNSIGNED_LONG_LONG) )
        return LLONG;

   return -1;
}

/******************************************************************************/
/** get_size() returns the number of bytes in MPI_Datatype argument
 *  @param    MPI_Datatype  t
 *  @return   number of bytes in t
 */
int get_size (MPI_Datatype t)
{
    if ( ( t == MPI_BYTE) ||
         ( t == MPI_CHAR) ||
         ( t == MPI_SIGNED_CHAR) ||
         ( t == MPI_UNSIGNED_CHAR) )
        return sizeof(char);

    if ( t == MPI_SHORT )
        return sizeof(short int);

    if ( t == MPI_DOUBLE )
        return sizeof(double);
    if (t == MPI_FLOAT)
        return sizeof(float);

    if ( ( t == MPI_INT) ||
         ( t == MPI_UNSIGNED) )
        return sizeof(int);

    if ( ( t == MPI_LONG) ||
         ( t == MPI_UNSIGNED_LONG) )
        return sizeof(long int);

    if ( ( t == MPI_LONG_LONG) ||
         ( t == MPI_LONG_LONG_INT) ||
         ( t == MPI_UNSIGNED_LONG_LONG) )
        return sizeof(long long int);

   return -1;
}