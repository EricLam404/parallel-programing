# Parallel Programming

This repository is dedicated to the CSCI 49365, Parallel Programming course at Hunter College, focusing on MPI (Message Passing Interface) assignments.

## Assignments

1. **Assignment 1: estimate_ln.c**
   - Description: Estimates the natural log of a number based on the number of segments
   - Build: mpicc -Wall -g -o estimate_ln estimate_ln.c -lm
   - Run: mpirun -H [hosts] estimate_ln [input_value] [num_segments]


