/**
 * @file    parallel_sort.cpp
 * @author  Patrick Flick <patrick.flick@gmail.com>
 * @brief   Implements the parallel, distributed sorting function.
 *
 * Copyright (c) 2014 Georgia Institute of Technology. All Rights Reserved.
 */

#include "parallel_sort.h"

// implementation of your parallel sorting
void parallel_sort(int * begin, int* end, MPI_Comm comm) {
    int  size, rank, err;

    err = MPI_Comm_size(comm, &size); ERR(err);
    err = MPI_Comm_rank(comm, &rank); ERR(err);

    int local_size = end - begin;
    int k = generate_random_number();
    // Right now assuming everyone has n/p elements
    int is_pivot_holder = (k/local_size) == rank;
    int pivot;
    
    if(is_pivot_holder){
        pivot = *(begin + (k % local_size));
        err = MPI_Bcast(&pivot, 1, MPI_INT, rank, comm); ERR(err);    
    }
    else{
        err = MPI_Bcast(&pivot, 1, MPI_INT, rank, comm); ERR(err);
    }

    //count elements lower and higher than pivot
    int local_low=0, local_high=0;
    for (int i = 0; i < local_size; i++)
    {
        if(begin[i] <= pivot)
            local_low++;
        else
            local_high++;
    }
    
    //store lower elements in local_lows and higher in local_highs
    int *local_lows = (int *)malloc(sizeof(int) * local_low);
    int *local_highs = (int *)malloc(sizeof(int) * local_high);
    int l = h = 0;
    for (int i = 0; i < local_size; i++){
        if(begin[i] <= pivot)
            local_lows[l++];
        else
            local_highs[h++];
    }

    //Gather total lows and total highs
    int *all_lows = (int *)malloc(sizeof(int) * size);
    MPI_Allgather(&low, 1, MPI_INT, all_lows, 1, MPI_INT, comm);

    int *all_highs = (int *)malloc(sizeof(int) * size);
    MPI_Allgather(&high, 1, MPI_INT, all_highs, 1, MPI_INT, comm);

    int total_low = 0;
    for (int i = 0; i < size; i++)
        total_low += all_lows[i];
    
    int total_high = 0;
    for (int i = 0; i < size; i++)
        total_high += all_highs[i];

    //Compute the communicator split based on total_highs, total_lows
    int lsize = (total_low * size)/ (total_high + total_low);
    if(!lsize)
        lsize++;
    int hsize = size - lsize;

    
    //coloring to split the communicator
    int color = (rank >= lsize);
    //prepare for alltoall transfer
    if(!color){
		//I will send local_highs to all processors with color 1
		int *sendcounts = (int *)malloc(sizeof(int) * size);
		int *sdispls = (int *)malloc(sizeof(int) * size);
		int *recvcounts = (int *)malloc(sizeof(int) * size);
		int *rdispls = (int *)malloc(sizeof(int) * size);
		for (int i = 0; i < size; i++){
			if(i < lsize)
				sendcounts[i] = 0;
			else
				sendcounts[i] = local_highs / hsize;
		}
		sendcounts[size - 1] += local_highs % hsize;
		sdispls[0] = 0;
		for (int i = 1; i < size; i++)
			sdispls[i] = sdispls[i - 1] + sendcounts[i];
		
		int *send_buf = (int *)malloc(sizeof(int) *)
	}
    //int new_rank = color ? (rank - (size / 2)) : (rank + (size / 2));
    

}

int generate_random_number(int m){
    srand48(43);
    long rand = lrand48();
    return rand % m;
}
/*********************************************************************
 *             Implement your own helper functions here:             *
 *********************************************************************/

// ...
