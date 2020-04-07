#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <bits/stdc++.h>
using namespace std;
/*
Debug code runs just if the flag is on; Setting it to 0 by default
*/
#define debug 0

//Generates random numbers array
void populate_random_array(int N, long int C, int rank, int size, double* arr){
    long int seed = C + rank;
    srand48(seed);
    int n = N / size;
    for (int i = 0; i < n; i++){
       if(debug)
            arr[i] = (rank*n + i+1);
        else
            arr[i] = drand48();
    }
    return;
}

int main(int argc, char** argv) {
  
    MPI_Init(&argc, &argv);

    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int N, C;
    //braodcasting N,C
    int payload[2];
    
    /*Assuming rank 0 is parsing the inputs and broadcasting*/
    if(!rank){
        if (argc != 3)
        {
            printf("Please provide the values of N and C");
            MPI_Finalize();
            return 1;
        }
        else{
            N = atoi(argv[1]);
            C = atoi(argv[2]);
            payload[0] = N;
            payload[1] = C;
        }
        if(debug)
            cout << "rank " << rank << " "<< N << " " << C << endl;
    }
    
    //Braodcasting the N,C values to all the processors
    MPI_Bcast(payload, 2, MPI_INT, 0, MPI_COMM_WORLD);

    double time;
    //Start the time clock
    time = MPI_Wtime();
    
    //Populate a random array 
    N = payload[0];
    C = payload[1];
    int sz = N/size;
    double *arr;
    arr = (double *)calloc(sz, sizeof(double));
    populate_random_array(N, C, rank, size, arr);
    
    //Computing local sum of the random array
    double local_sum = 0.0;
    for(int i=0; i<sz; i++)
        local_sum += arr[i];
    
    //Algorithm as discussed in the class
    int j = 1;
    while ( j < size)
    {
    int dst = rank ^ j;
        if((rank & j) != 0){
            if(debug)
               cout << "Src " << rank << " Des "<< dst << endl;
            MPI_Send(&local_sum, 1, MPI_DOUBLE, dst, 99, MPI_COMM_WORLD);
            break;
        }
        else{
            double tmp;
            MPI_Status status;
            MPI_Recv(&tmp, 1, MPI_DOUBLE, dst, 99, MPI_COMM_WORLD, &status);
            local_sum += tmp;
        }
        j = j << 1;
    }
    
    //free the array, as we no longer need it
    free(arr);
    
    time = MPI_Wtime() - time;
    //using MPI_reduce to calculate max of all the times
    double max_time;
    MPI_Reduce(&time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        
    if(!rank){
        printf("N= %d, P= %d, C= %d, S= %lf \nTime= %12.5e\n", N, size, C, local_sum, max_time);
    }
    MPI_Finalize();
    return 0;
}
