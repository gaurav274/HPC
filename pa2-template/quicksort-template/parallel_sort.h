/**
 * @file    parallel_sort.h
 * @author  Patrick Flick <patrick.flick@gmail.com>
 * @brief   Declares the parallel sorting function.
 *
 * Copyright (c) 2014 Georgia Institute of Technology. All Rights Reserved.
 */

#ifndef PARALLEL_SORT_H
#define PARALLEL_SORT_H

#define ERR(err) if (err) return err
#define DEBUG(rank, msg) std::cout << rank << ":" << msg << std::endl;
#include <mpi.h>
#include <stdlib.h>
#include <iostream> 
#include <cstring>
#include <algorithm>
#include <unistd.h>
using namespace std;
/**
 * @brief   Parallel, distributed sorting over all processors in `comm`. Each
 *          processor has the local input [begin, end).
 *
 * Note that `end` is given one element beyond the input. This corresponds to
 * the API of C++ std::sort! You can get the size of the local input with:
 * int local_size = end - begin;
 *
 * @param begin Pointer to the first element in the input sequence.
 * @param end   Pointer to one element past the input sequence. Don't access this!
 * @param comm  The MPI communicator with the processors participating in the
 *              sorting.
 */
void parallel_sort(int * begin, int* end, MPI_Comm comm);


/**
 * @brief   Helper function to recursively sort all the number in the comm
 *
 *
 * @param local_arr     Pointer to the first element in the local sequence.
 * @param local_size    Size of the local_arr
 * @param sorted_local_arr     Pointer to the first element in the resulting sorted local sequence.
 * @param sorted_local_size    Size of the sorted_local_arr
 * @param comm  The MPI communicator with the processors participating in the
 *              sorting.
 */
int parallel_sort_recursive(int *local_arr, int local_size, int **sorted_local_arr, int *sorted_local_size, MPI_Comm comm);


//Random number generator
int generate_random_number(int m);

//Generates Pivot
int getPivot(int *local_arr, int local_size, MPI_Comm comm);

/*********************************************************************
 *              Declare your own helper functions here               *
 *********************************************************************/

// ...

#endif // PARALLEL_SORT_H
