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

/*
 * TODO: Implement your solutions here
 */

void distribute_vector(const int n, double *input_vector, double **local_vector, MPI_Comm comm)
{
    if (get_coord_by_dim(comm, 1) != 0)
        return;

    int belongs[2] = {1, 0};
    MPI_Comm commcol;
    MPI_Cart_sub(comm, belongs, &commcol);
    int col_rank, col_size;
    MPI_Comm_rank(commcol, &col_rank);
    MPI_Comm_size(commcol, &col_size);
    int local_x_size = block_decompose_by_dim(n, comm, 0);
    int sendcounts[col_size], displs[col_size];
    if (!col_rank)
    {
        for (int i = 0; i < col_size; i++)
            sendcounts[i] = block_decompose(n, col_size, i);
        displs[0] = 0;
        for (int i = 1; i < col_size; i++)
            displs[i] = displs[i - 1] + sendcounts[i - 1];
    }
    MPI_Scatterv(input_vector, sendcounts, displs, MPI_DOUBLE, *local_vector, local_x_size, MPI_DOUBLE, 0, commcol);

    MPI_Comm_free(commcol);
}

// gather the local vector distributed among (i,0) to the processor (0,0)
void gather_vector(const int n, double *local_vector, double *output_vector, MPI_Comm comm)
{
    if (get_coord_by_dim(comm, 1) != 0)
        return;

    int belongs[2] = {1, 0};
    MPI_Comm commcol;
    MPI_Cart_sub(comm, belongs, &commcol);
    int col_rank, col_size;
    MPI_Comm_rank(commcol, &col_rank);
    MPI_Comm_size(commcol, &col_size);
    int local_x_size = block_decompose_by_dim(n, comm, 0);
    int recvcounts[col_size], displs[col_size];
    if (!col_rank)
    {
        for (int i = 0; i < col_size; i++)
            recvcounts[i] = block_decompose(n, col_size, i);
        displs[0] = 0;
        for (int i = 1; i < col_size; i++)
            displs[i] = displs[i - 1] + sendcounts[i - 1];
    }
    MPI_Gatherv(local_vector, local_x_size, MPI_DOUBLE, output_vector, recvcounts, displs, MPI_DOUBLE, 0, commcol);

    MPI_Comm_free(commcol);

}

void distribute_matrix(const int n, double *input_matrix, double **local_matrix, MPI_Comm comm)
{
    int belongs[2] = {1, 0};
    MPI_Comm commcol;
    MPI_Cart_sub(comm, belongs, &commcol);
    int col_rank, col_size;
    MPI_Comm_rank(commcol, &col_rank);
    MPI_Comm_size(commcol, &col_size);
    int local_x_size = n*block_decompose_by_dim(n, comm, 0);
    double *local_vector = new double[local_x_size];
    int sendcounts[col_size], displs[col_size];
    if(get_coord_by_dim(comm, 1) == 0){
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

    //scatter from col=0
    int belongs[2] = {0, 1};
    MPI_Comm commrow;
    MPI_Cart_sub(comm, belongs, &commrow);
    int row_rank, row_size;
    MPI_Comm_rank(commcol, &row_rank);
    MPI_Comm_size(commcol, &row_size);
    int local_x_size = block_decompose_by_dim(n, comm, 1);
    int local_y_size = block_decompose_by_dim(n, comm, 0);
    //local_vector = new double[local_x_size];
    int sendcounts[row_size], displs[row_size];
    if (!row_rank){
        for (int i = 0; i < row_size; i++)
            sendcounts[i] = block_decompose(n, row_size, i);
        displs[0] = 0;
        for (int i = 1; i < row_size; i++)
            displs[i] = displs[i - 1] + sendcounts[i - 1];
        }
    for (int i = 0; i < local_y_size; i++)
        MPI_Scatterv(local_vector + i*n, sendcounts, displs, MPI_DOUBLE, *local_matrix + local_x_size*i, local_x_size, MPI_DOUBLE, 0, commrow);

    MPI_Comm_free(commcol);
    MPI_Comm_free(commrow);
    free(local_vector);
}

void transpose_bcast_vector(const int n, double *col_vector, double *row_vector, MPI_Comm comm)
{
    int dims[2];
    int periods[2];
    int coords[2];
    MPI_Cart_get(comm, 2, dims, periods, coords);
    int col = coord[0];
    int local_size = block_decompose_by_dim(n, comm, 0);
    /* Create 1D row subgrids */
    int belongs[2];
    MPI_Comm commrow;
    belongs[0] = 0;
    belongs[1] = 1; // this dimension belongs to subgrid
    MPI_Cart_sub(comm, belongs, &commrow);
    int row_rank;
    int rank;
    MPI_Comm_rank(commrow, &row_rank);
    if (!row_rank)
    {
        MPI_Send(col_vector, local_size, MPI_DOUBLE, col, '99', commrow);
    }
    if (row_rank == col)
    {
        MPI_Send(row_vector, local_size, MPI_DOUBLE, 0, '99', commrow, MPI_STATUS_IGNORE);
    }

    // MPI_Barrier(comm);

    belongs[0] = 1;
    belongs[1] = 0;
    MPI_Cart_sub(comm, belongs, &commcol);
    MPI_Bcast(row_vector, local_size, MPI_DOUBLE, col, commcol);

    MPI_Comm_free(&commrow);
    MPI_Comm_free(&commcol);
}

void distributed_matrix_vector_mult(const int n, double *local_A, double *local_x, double *local_y, MPI_Comm comm)
{

    int local_x_size = block_decompose_by_dim(n, comm, 0);
    int local_y_size = block_decompose_by_dim(n, comm, 1);
    transpose_bcast_vector(n, local_x, new_local_x, comm);

    for (int row = 0; row < local_y_size; row++)
    {
        for (int col = 0; col < local_x_size; col++)
            local_y[row] += local_A[row * local_x_size + col] * new_local_x[col];
    }

    int belongs[2];
    MPI_Comm commrow;
    belongs[0] = 0;
    belongs[1] = 1; // this dimension belongs to subgrid
    MPI_Cart_sub(comm, belongs, &commrow);
    int row_rank;
    MPI_Comm_rank(commrow, &row_rank);
    MPI_Reduce(local_y, local_y, local_y_size, MPI_SUM, 0, commrow);
    MPI_Comm_free(&commrow);
}

// Solves Ax = b using the iterative jacobi method
void distributed_jacobi(const int n, double *local_A, double *local_b, double *local_x,
                        MPI_Comm comm, int max_iter, double l2_termination)
{
    // TODO
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
