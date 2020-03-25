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
	//while(1);
	int size, rank, err;
	int *new_arr = NULL;
	int new_size, old_size = end - begin;
	
	err = MPI_Comm_size(comm, &size); 
    err = MPI_Comm_rank(comm, &rank); 

	parallel_sort_recursive(begin, old_size, &new_arr, &new_size, comm);
   
    //restored based on original size

    //prepare for alltoall transfer
    int *sendcounts = (int *)malloc(sizeof(int) * size);
	int *sdispls = (int *)malloc(sizeof(int) * (size+1));
	int *recvcounts = (int *)malloc(sizeof(int) * size);
	int *rdispls = (int *)malloc(sizeof(int) * (size+1));
	//int *send_buf, *recv_buf;

	int *olddispls = (int *)malloc(sizeof(int) * (size+1));
	int *newdispls = (int *)malloc(sizeof(int) * (size+1));

	int *old_sizes = (int *)malloc(sizeof(int) * size);
	err = MPI_Allgather(&old_size, 1, MPI_INT, old_sizes, 1, MPI_INT, comm); //ERR(err);

	int *new_sizes = (int *)malloc(sizeof(int) * size);
	err = MPI_Allgather(&new_size, 1, MPI_INT, new_sizes, 1, MPI_INT, comm); //ERR(err);

	olddispls[0] = 0;
	for (int i = 1; i <= size; i++)
		olddispls[i] = olddispls[i - 1] + old_sizes[i-1];

	
	newdispls[0] = 0;
	for (int i = 1; i <= size; i++)
		newdispls[i] = newdispls[i - 1] + new_sizes[i-1];

	
	for (int i = 0; i < size; i++) {
		if (newdispls[i+1] <= olddispls[rank] || olddispls[rank+1] <= newdispls[i]) {
			recvcounts[i] = 0;
		}
		else{
			recvcounts[i] = std::min(olddispls[rank+1],newdispls[i+1]) - std::max(olddispls[rank],newdispls[i]);
		}
	}

	for (int i = 0; i < size; i++) {
		if (olddispls[i+1] <= newdispls[rank] || newdispls[rank+1] <= olddispls[i]) {
			sendcounts[i] = 0;
		}
		else{
			sendcounts[i] = std::min(newdispls[rank+1],olddispls[i+1]) - std::max(newdispls[rank],olddispls[i]);
		}
	}

	sdispls[0] = 0;
	for (int i = 1; i <= size; i++)
		sdispls[i] = sdispls[i - 1] + sendcounts[i-1];

	
	rdispls[0] = 0;
	for (int i = 1; i <= size; i++)
		rdispls[i] = rdispls[i - 1] + recvcounts[i-1];

	err = MPI_Alltoallv(new_arr, sendcounts, sdispls, MPI_INT, begin, recvcounts, rdispls, MPI_INT, comm); //ERR(err);



}
int parallel_sort_recursive(int *local_arr, int local_size, int **sorted_local_arr, int *sorted_local_size, MPI_Comm comm){
	int  size, rank, err;

    err = MPI_Comm_size(comm, &size); ERR(err);
    err = MPI_Comm_rank(comm, &rank); ERR(err);
	
	if (size == 1)
	{
		sort(local_arr, local_arr + local_size);
		*sorted_local_size = local_size;
		*sorted_local_arr = local_arr;
		return 0;
	}
	//int k = generate_random_number(size);
	// Right now assuming everyone has n/p elements
    int pivot = getPivot(local_arr,local_size,comm);
    
    //count elements lower and higher than pivot
    int local_low=0, local_high=0;
    for (int i = 0; i < local_size; i++)
    {
        if(local_arr[i] <= pivot)
            local_low++;
        else
            local_high++;
    }
    
    //store lower elements in local_lows and higher in local_highs
    int *local_lows = (int *)malloc(sizeof(int) * local_low);
    int *local_highs = (int *)malloc(sizeof(int) * local_high);
    int l = 0,h = 0;
    for (int i = 0; i < local_size; i++){
        if(local_arr[i] <= pivot)
            local_lows[l++] = local_arr[i];
        else
            local_highs[h++] = local_arr[i];
    }
    
    
    //Gather total lows and total highs
    int *all_lows = (int *)malloc(sizeof(int) * size);
    err = MPI_Allgather(&local_low, 1, MPI_INT, all_lows, 1, MPI_INT, comm); ERR(err);

    int *all_highs = (int *)malloc(sizeof(int) * size);
    err = MPI_Allgather(&local_high, 1, MPI_INT, all_highs, 1, MPI_INT, comm); ERR(err);

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
    if(!hsize){
    	hsize++;
    	lsize--;
    }
    //coloring to split the communicator
    int color = (rank >= lsize);

    //prepare for alltoall transfer
    int *sendcounts = (int *)malloc(sizeof(int) * size);
	int *sdispls = (int *)malloc(sizeof(int) * size);
	int *recvcounts = (int *)malloc(sizeof(int) * size);
	int *rdispls = (int *)malloc(sizeof(int) * size);
	int *send_buf, *recv_buf;
	if (!color)
	{
		//I will send local_highs to all processors with color 1
		//and recieve local_lows from all processors with color 1
		
		//Setting up sending logic
		for (int i = 0; i < size; i++){
			if(i < lsize)
				sendcounts[i] = 0;
			else
				sendcounts[i] = local_high / hsize;
		}
		sendcounts[size - 1] += local_high % hsize;
		
		send_buf = local_highs;
		//Setting up receiving logic
		for (int i = 0; i < size; i++){
			if(i < lsize)
				recvcounts[i] = 0;
			else
				recvcounts[i] = all_lows[i] / lsize + ((rank==lsize-1) ? all_lows[i] % lsize : 0);
			
		}		
	}
	else{
		//I will send local_lows to all processors with color 0
		//recieve local_highs from them

		//Setting up sending logic
		for (int i = 0; i < size; i++){
			if(i >= lsize)
				sendcounts[i] = 0;
			else
				sendcounts[i] = local_low / lsize;
		}
		sendcounts[lsize - 1] += local_low % lsize;
		
		send_buf = local_lows;
		
		//Setting up receiving logic
		for (int i = 0; i < size; i++){
			if(i >= lsize)
                recvcounts[i] = 0;
			else
				recvcounts[i] = all_highs[i] / hsize + ((rank==size-1) ? all_highs[i] % hsize : 0);
		}
	}
	
	sdispls[0] = 0;
	for (int i = 1; i < size; i++)
		sdispls[i] = sdispls[i - 1] + sendcounts[i-1];

	
	rdispls[0] = 0;
	for (int i = 1; i < size; i++)
		rdispls[i] = rdispls[i - 1] + recvcounts[i-1];
    
    //create recv_buf to recv elements from other processors and local elements
    int recv_buf_size; 
	if(!color)
        recv_buf_size = rdispls[size - 1] + recvcounts[size - 1] + local_low;
    else
        recv_buf_size = rdispls[size - 1] + recvcounts[size - 1] + local_high;
    
    recv_buf = (int *)malloc(sizeof(int) * recv_buf_size);

    err = MPI_Alltoallv(send_buf, sendcounts, sdispls, MPI_INT, recv_buf, recvcounts, rdispls, MPI_INT, comm);
	ERR(err);
    
    //Append local elements to the recv_buf
    
    if(!color){
        for(int i=0;i<local_low;i++)
            recv_buf[rdispls[size - 1] + recvcounts[size - 1] + i] = local_lows[i];
           
    }
    else{
         for(int i=0;i<local_high;i++)
            recv_buf[rdispls[size - 1] + recvcounts[size - 1] + i] = local_highs[i];
    }
   
    MPI_Comm subcomm;
	err = MPI_Comm_split(comm, color, rank, &subcomm); ERR(err);
	
	int *new_local_arr = NULL;
	int new_local_size;
    
	err = parallel_sort_recursive(recv_buf, recv_buf_size, &new_local_arr, &new_local_size, subcomm);
	ERR(err);

	int *tmp = (int *)malloc(sizeof(int) * new_local_size);
	memcpy(tmp, new_local_arr, sizeof(int) * new_local_size);
	*sorted_local_size = new_local_size;
	*sorted_local_arr = tmp;
	return 0;
}
int generate_random_number(int m){
    srand48(43);
    long rand = lrand48();
    return rand % m;
}
int getPivot(int *local_arr, int local_size, MPI_Comm comm){
	int  size, rank, err;

    err = MPI_Comm_size(comm, &size); ERR(err);
    err = MPI_Comm_rank(comm, &rank); ERR(err);

    int pivot = 0;
    sort(local_arr, local_arr + local_size);
    if(local_size)	pivot = local_arr[local_size/2];
    
    int *lmedians = (int *)malloc(sizeof(int) * size);
    err = MPI_Allgather(&pivot, 1, MPI_INT, lmedians, 1, MPI_INT, comm); ERR(err);

    sort(lmedians, lmedians + size);
    return lmedians[size/2];
}
/*********************************************************************
 *             Implement your own helper functions here:             *
 *********************************************************************/

// ...
