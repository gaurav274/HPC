/**
 * @file    mpi_jacobi.cpp
 * @author  Patrick Flick <patrick.flick@gmail.com>
 * @brief   Implements MPI functions for distributing vectors and matrixes,
 *          parallel distributed matrix-vector multiplication and Jacobi's
 *          method.
 *
 * Copyright (c) 2014 Georgia Institute of Technology. All Rights Reserved.
 */

#include "mpi_jacobi.h"
#include "jacobi.h"
#include "utils.h"

#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <vector>
#include <string.h>

/*
 * TODO: Implement your solutions here
 */

void distribute_vector(const int n, double *input_vector, double **local_vector, MPI_Comm comm)
{
    int rank;
    MPI_Comm_rank(comm, &rank);
    //Get column communicator
    int belongs[2] = {1, 0};
    MPI_Comm commcol;
    MPI_Cart_sub(comm, belongs, &commcol);
    int col_rank, col_size;
    MPI_Comm_rank(commcol, &col_rank);
    MPI_Comm_size(commcol, &col_size);
    
    MPI_Comm_rank(comm, &rank); 
    // only (i,0) participates
    if (get_coord_by_dim(comm, 1) != 0)
        return;

    //Get the local size
    int local_x_size = block_decompose_by_dim(n, comm, 0);
    // Scatter from the rank 0 processor
    int sendcounts[col_size], displs[col_size];
    for (int i = 0; i < col_size; i++)
        sendcounts[i] = block_decompose(n, col_size, i);
    displs[0] = 0;
    for (int i = 1; i < col_size; i++)
        displs[i] = displs[i - 1] + sendcounts[i - 1];
    *local_vector = new double[local_x_size];
    MPI_Scatterv(input_vector, sendcounts, displs, MPI_DOUBLE, *local_vector, local_x_size, MPI_DOUBLE, 0, commcol);
    
    MPI_Comm_free(&commcol);
}

// gather the local vector distributed among (i,0) to the processor (0,0)
void gather_vector(const int n, double *local_vector, double *output_vector, MPI_Comm comm)
{
    //Get column communicator
    int belongs[2] = {1, 0};
    MPI_Comm commcol;
    MPI_Cart_sub(comm, belongs, &commcol);
    int col_rank, col_size;
    MPI_Comm_rank(commcol, &col_rank);
    MPI_Comm_size(commcol, &col_size);
    
    int local_x_size = block_decompose_by_dim(n, comm, 0);
    int recvcounts[col_size], displs[col_size];
    
    // only (i,0) participates
    if (get_coord_by_dim(comm, 1) != 0)
        return;

    for (int i = 0; i < col_size; i++)
        recvcounts[i] = block_decompose(n, col_size, i);
    displs[0] = 0;
    for (int i = 1; i < col_size; i++)
        displs[i] = displs[i - 1] + recvcounts[i - 1];
    MPI_Gatherv(local_vector, local_x_size, MPI_DOUBLE, output_vector, recvcounts, displs, MPI_DOUBLE, 0, commcol);

    MPI_Comm_free(&commcol);

}

void distribute_matrix(const int n, double *input_matrix, double **local_matrix, MPI_Comm comm)
{
    //Get column communicator
    int belongs[2] = {1, 0};
    MPI_Comm commcol;
    MPI_Cart_sub(comm, belongs, &commcol);
    int col_rank, col_size;
    MPI_Comm_rank(commcol, &col_rank);
    MPI_Comm_size(commcol, &col_size);
    
    int local_x_size = n*block_decompose_by_dim(n, comm, 0);
    double *local_vector = new double[local_x_size];
    
    //distributing matrix along first column
    if(get_coord_by_dim(comm, 1) == 0){
        int sendcounts[col_size], displs[col_size];
        if (!col_rank)
        {
            for (int i = 0; i < col_size; i++)
                sendcounts[i] = n*block_decompose(n, col_size, i);
            displs[0] = 0;
            for (int i = 1; i < col_size; i++)
                displs[i] = displs[i - 1] + sendcounts[i - 1];
        }
        MPI_Scatterv(input_matrix, sendcounts, displs, MPI_DOUBLE, local_vector, local_x_size, MPI_DOUBLE, 0, commcol);
    }
    
    //Get row communicator
    belongs[0] = 0;
    belongs[1] = 1;
    MPI_Comm commrow;
    MPI_Cart_sub(comm, belongs, &commrow);
    int row_rank, row_size;
    MPI_Comm_rank(commrow, &row_rank);
    MPI_Comm_size(commrow, &row_size);
    
    local_x_size = block_decompose_by_dim(n, comm, 1);
    int local_y_size = block_decompose_by_dim(n, comm, 0);
    
    *local_matrix = new double[local_x_size* local_y_size];
    int sendcounts[row_size], displs[row_size];
    
    //distributing matrix along every row
    if (!row_rank){
        for (int i = 0; i < row_size; i++)
            sendcounts[i] = block_decompose(n, row_size, i);
        displs[0] = 0;
        for (int i = 1; i < row_size; i++)
            displs[i] = displs[i - 1] + sendcounts[i - 1];
        }
    for (int i = 0; i < local_y_size; i++)
        MPI_Scatterv(local_vector + i*n, sendcounts, displs, MPI_DOUBLE, *local_matrix + local_x_size*i, local_x_size, MPI_DOUBLE, 0, commrow);
    
    int rank;
    MPI_Comm_rank(comm, &rank); 
    MPI_Comm_free(&commcol);
    MPI_Comm_free(&commrow);
    free(local_vector);
}

void transpose_bcast_vector(const int n, double *col_vector, double *row_vector, MPI_Comm comm)
{
    int row = get_coord_by_dim(comm, 0);
    int local_size = block_decompose_by_dim(n, comm, 0);
    /* Create 1D row subgrids */
    int belongs[2] = {0,1};
    MPI_Comm commrow;
    MPI_Cart_sub(comm, belongs, &commrow);
    int row_rank;
    int rank;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_rank(commrow, &row_rank);

    // nothing to do for (0,0)
    if(row != 0){
        if (!row_rank)
            MPI_Send(col_vector, local_size, MPI_DOUBLE, row, 99, commrow);
        if (row_rank == row)
            MPI_Recv(row_vector, local_size, MPI_DOUBLE, 0, 99, commrow, MPI_STATUS_IGNORE);
    }
    if(!row && !row_rank)
        memcpy(row_vector, col_vector, local_size*sizeof(double));

    belongs[0] = 1;
    belongs[1] = 0;
    MPI_Comm commcol;
    int local_size_y = block_decompose_by_dim(n, comm, 1);
    int col = get_coord_by_dim(comm, 1);
    MPI_Cart_sub(comm, belongs, &commcol);
    MPI_Bcast(row_vector, local_size_y, MPI_DOUBLE, col, commcol);

    MPI_Comm_free(&commrow);
    MPI_Comm_free(&commcol);
}

void distributed_matrix_vector_mult(const int n, double *local_A, double *local_x, double *local_y, MPI_Comm comm)
{
    int rank;
    MPI_Comm_rank(comm, &rank); 
    
    int local_y_size = block_decompose_by_dim(n, comm, 0);
    int local_x_size = block_decompose_by_dim(n, comm, 1);
    double *new_local_x = new double[local_x_size];
    transpose_bcast_vector(n, local_x, new_local_x, comm);
    
    double *new_local_y = new double[local_y_size];
    
    //local computation for y = Ax
    for (int row = 0; row < local_y_size; row++)
    {
        new_local_y[row] = 0.0;
        for (int col = 0; col < local_x_size; col++)
            new_local_y[row] += local_A[row * local_x_size + col] * new_local_x[col];
    }

    int belongs[2];
    MPI_Comm commrow;
    belongs[0] = 0;
    belongs[1] = 1; // this dimension belongs to subgrid
    MPI_Cart_sub(comm, belongs, &commrow);
    int row_rank;
    MPI_Comm_rank(commrow, &row_rank);
    
    //reducing result along the row which will be owned by processor in first column
    MPI_Reduce(new_local_y, local_y, local_y_size, MPI_DOUBLE, MPI_SUM, 0, commrow);
    
    MPI_Comm_free(&commrow);
    free(new_local_x);
    free(new_local_y);
}

// Solves Ax = b using the iterative jacobi method
void distributed_jacobi(const int n, double *local_A, double *local_b, double *local_x,
                        MPI_Comm comm, int max_iter, double l2_termination)
{
    int rank;
    MPI_Comm_rank(comm, &rank);
    int local_x_size = block_decompose_by_dim(n, comm, 0);
    int local_y_size = block_decompose_by_dim(n, comm, 1);
    int row = get_coord_by_dim(comm, 0);
    int col = get_coord_by_dim(comm, 1);
    double diagnol[local_y_size];
    double local_diagnol[local_y_size];
    double local_R[local_y_size*local_x_size];
    memcpy(local_R, local_A, sizeof(double)*local_y_size*local_x_size);
    
    // only row = col will have diagnol elements
    // construct R and make local_diagnol
    int d=0;
    if(row == col){
        for(int i=0;i<local_x_size; i++){
            for(int j=0;j<local_y_size; j++){
                if(i == j){
                    local_diagnol[d++] = local_A[i*local_x_size+ j];
                    local_R[i*local_x_size+ j] = 0.0;
                }
            }
        }
    }
    
    int belongs[2] = {0,1};
    MPI_Comm commrow;
    MPI_Cart_sub(comm, belongs, &commrow);
    belongs[0] = 1; belongs[1] = 0;
    MPI_Comm commcol;
    MPI_Cart_sub(comm, belongs, &commcol);
    
    int row_rank;
    MPI_Comm_rank(commrow, &row_rank);
    
    // rank(0,0) copies locally
    memcpy(diagnol, local_diagnol, local_y_size*sizeof(double));
    // others communicates
    if(row != 0){
        // Send local diagnol to first column
        if(row == col){
            MPI_Send(local_diagnol, local_y_size, MPI_DOUBLE, 0, 99, commrow);
        }
        //first col recieves
        if(col == 0){
            MPI_Recv(diagnol, local_y_size, MPI_DOUBLE, row, 99, commrow, MPI_STATUS_IGNORE);
        }
    }
    
    // set local_x = 0 
    for(int i=0;i<local_y_size;i++)
        local_x[i] = 0.0;
    
    // Iteration
    int terminate = 0;
    for(int iter=0;iter<max_iter && !terminate;iter++){
        double *local_y = new double[local_y_size];
        distributed_matrix_vector_mult(n, local_R, local_x, local_y, comm);

        //only do in first columns
        if(!col){
            for(int i=0; i< local_y_size; i++)
                local_x[i] = (local_b[i] - local_y[i])/diagnol[i];
        }
        // ||A*x - b||
        distributed_matrix_vector_mult(n, local_A, local_x, local_y, comm);
        double local_norm=0, global_norm;
        if(!col){
            for(int i=0; i< local_y_size; i++)
                local_norm += (local_b[i] - local_y[i])*(local_b[i] - local_y[i]);
        }
        
        MPI_Reduce(&local_norm, &global_norm, 1, MPI_DOUBLE, MPI_SUM, 0, commcol);
        MPI_Bcast(&global_norm, 1, MPI_DOUBLE, 0, comm);

        if(sqrt(global_norm) <= l2_termination){
            // we need to terminate so broadcast everyone to stop
            terminate = 1;
        }
        free(local_y);
    }
    //MPI_Send(local_diagnol, local_y_size, MPI_DOUBLE, 0, 99, commrow);
}

// wraps the distributed matrix vector multiplication
void mpi_matrix_vector_mult(const int n, double *A,
                            double *x, double *y, MPI_Comm comm)
{
    // distribute the array onto local processors!
    double *local_A = NULL;
    double *local_x = NULL;
    
    distribute_matrix(n, &A[0], &local_A, comm);
    distribute_vector(n, &x[0], &local_x, comm);
    
    
    // allocate local result space
    double *local_y = new double[block_decompose_by_dim(n, comm, 0)];
    distributed_matrix_vector_mult(n, local_A, local_x, local_y, comm);

    // gather results back to rank 0
    gather_vector(n, local_y, y, comm);
}

// wraps the distributed jacobi function
void mpi_jacobi(const int n, double *A, double *b, double *x, MPI_Comm comm,
                int max_iter, double l2_termination)
{
    // distribute the array onto local processors!
    double *local_A = NULL;
    double *local_b = NULL;
    distribute_matrix(n, &A[0], &local_A, comm);
    distribute_vector(n, &b[0], &local_b, comm);

    // allocate local result space
    double *local_x = new double[block_decompose_by_dim(n, comm, 0)];
    distributed_jacobi(n, local_A, local_b, local_x, comm, max_iter, l2_termination);

    // gather results back to rank 0
    gather_vector(n, local_x, x, comm);
}
