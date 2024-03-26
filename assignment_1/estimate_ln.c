/******************************************************************************
Title : estimate_ln.c
Author : Eric Lam
Created on : March 2, 2024
Description : Estimates the natural log of a number based on the number of 
              segments
Purpose : Demonstrates parallel programing using MPI by splitting the estimates
of each local sum to a processor and then adding it into a total sum
Usage : estimates natural log of a number
Build with : mpicc -Wall -g -o estimate_ln estimate_ln.c -lm
Modifications: March 6, 2024
Added Usage Error handling
******************************************************************************/

#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include <stdbool.h>
#include <math.h>
#include <stdlib.h>
#include "mpi.h"

#define ROOT 0

/** is_number()
 * checks if the character array is a number
 * returns true if it is a number and false if not 
*/
bool is_number(char number[]);

/** approximate_ln()
 *  This returns an approximation of the natural log of an input value
 *  To approximate ln we can use
 *  d/dx(ln(x)) = 1/x 
 *  and compute the area under the curve of 1/x in the interval from 1 to x
 *  We use the midpoint reimann sum to approximate the area under the curve. 
 *  We divide [1, x] into n segments of length (x+1)/n each, and compute the 
 *  value of 1/x at the midpoint of each segment. By summing these
 *  values and then multiplying by the length of each segment, we get
 *  the total area under the curve from [1, x] and an approximate of ln(x).
*/
double approximate_ln ( double input_value, int num_segments, int id, int p);

int main( int argc, char *argv[] )
{
    int id;                 /* rank of executing process    */
    int p;                  /* number of processes          */
    double ln_estimate;     /* estimated value of ln(x)     */
    double local_ln;        /* each process's contribution  */
    double input_value;     /* value to calculate the ln of */
    int num_intervals;      /* number of terms in series    */
    double elapsed_time;

    MPI_Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &id );
    MPI_Comm_size (MPI_COMM_WORLD, &p);

    MPI_Bcast(&input_value, 1, MPI_INT, ROOT, MPI_COMM_WORLD); 
    MPI_Bcast(&num_intervals, 1, MPI_INT, ROOT, MPI_COMM_WORLD); 
    MPI_Barrier (MPI_COMM_WORLD);
    elapsed_time = - MPI_Wtime();

    if(argc < 3){
        if(ROOT == id){
            printf("USAGE ERROR: insufficient command line arguments \n");
        }
        MPI_Finalize(); 
        return 0;
    }
    if(ROOT == 1){
        input_value = atof(argv[1]);
        num_intervals = atoi(argv[2]);
    }
    if (!is_number(argv[1]) || !is_number(argv[2])){
        if(ROOT == id){
            printf("USAGE ERROR: command line arguments are not numbers \n");
        }
    } else if (input_value < 1 || num_intervals < 1){
        if(ROOT == id){
            printf("USAGE ERROR: command line arguments are not numbers greater than 1 \n");
        }
    } else {
        local_ln = approximate_ln(input_value,num_intervals, id, p);
        MPI_Reduce(&local_ln, 
                &ln_estimate, 1, MPI_DOUBLE, 
                MPI_SUM, ROOT, MPI_COMM_WORLD);
        MPI_Barrier (MPI_COMM_WORLD);

        elapsed_time += MPI_Wtime();

        if (ROOT == id) {
            printf("%.16g    %.16f    %.16f    %8.6f seconds\n",
                input_value, ln_estimate, fabs(ln_estimate - log(input_value)), elapsed_time);
            fflush(stdout);
        }
    }
    MPI_Finalize(); 
    return 0;
}

bool is_number(char number[])
{
   /* initalizes dot_count to keep track of total decimal points */
    int dot_count = 0;

    /* moves the string one space to the right if it is negative */
    if (*number == '-') {
        number++;
    }
    /* loops through the whole string and checks to see if each 
       character is a digit and there are at most 1 deciaml point */
    while (*number) {
        if(*number == '.') {
            dot_count++;
            if (dot_count > 1) {
                return false;
            }
        } else if (!isdigit(*number)) {
            return false;
        }
        number++;
    }
    return true;
}

double approximate_ln ( double input_value, int num_segments, int id, int p)
{
    double dx, sum, x;
    int i;

    /* Set dx to the width of each segments */
    dx = (input_value - 1.0) / (double) num_segments;

    /* Initialize sum */
    sum = 0.0;

    /* Each process will compute its share of the segments. If the
       segments are numbered 1, 2, 3, ...,n, from left to right, then
       process id computes segment k if id = (k-1) % p, or equivalently
       it computes segments id+1, id+p+1, id+2p+1, ... up to id+mp+1, 
       where m is the largest number such that id+mp+1 <= num_segments. */
    for (i = id + 1; i <= num_segments; i += p) { 
        x = 1.0 + dx * ((double)i - 0.5);
        sum += 1.0 / x;  
    } 
    return  dx * sum; 
}