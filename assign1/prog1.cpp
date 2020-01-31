#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#define debug 0
#include <bits/stdc++.h>
using namespace std;
//Generates random numbers and compute their sum
double compute_local_sum(int N, long int C, int rank, int size){
    long int seed = C + rank;
    srand48(seed);
    int n = N / size;
    double sum = 0.0;
    for (int i = 0; i < n; i++){
        if(debug)
            sum += (rank*n + i+1);
        else
            sum += drand48();
    }
    return sum;
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

    MPI_Bcast(payload, 2, MPI_INT, 0, MPI_COMM_WORLD);

    double time;
    if(!rank)
	time = MPI_Wtime();
    double local_sum = compute_local_sum(payload[0], payload[1], rank, size);
    
    //if(debug)
       // cout << "rank " << rank << " "<< local_sum << endl;
 
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
    
    if(!rank){
        time = MPI_Wtime() - time;
        printf("N= %d, P= %d, C= %d, S= %lf \nTime= %12.5e\n", N, size, C, local_sum, time);
    }
    MPI_Finalize();
    return 0;
}
