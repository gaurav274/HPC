/**
 * @file    jacobi.cpp
 * @author  Patrick Flick <patrick.flick@gmail.com>
 * @brief   Implements matrix vector multiplication and Jacobi's method.
 *
 * Copyright (c) 2014 Georgia Institute of Technology. All Rights Reserved.
 */
#include "jacobi.h"

/*
 * TODO: Implement your solutions here
 */

// my implementation:
#include <iostream>
#include <math.h>
#include <cstring>
#include <stdlib.h>
#include <stdio.h>

// Calculates y = A*x for a square n-by-n matrix A, and n-dimensional vectors x
// and y
void matrix_vector_mult(const int n, const double *A, const double *x, double *y)
{
    for (int i = 0; i < n; i++)
    {
        double ls = 0;
        for (int j = 0; j < n; j++)
        {
            ls += (A[i*n + j] * x[j]);
        }
        y[i] = ls;
    }
}

// Calculates y = A*x for a n-by-m matrix A, a m-dimensional vector x
// and a n-dimensional vector y
void matrix_vector_mult(const int n, const int m, const double *A, const double *x, double *y)
{
    for (int i = 0; i < n; i++)
    {
        double ls = 0;
        for (int j = 0; j < m; j++)
        {
            ls += (A[i*m + j] * x[j]);
        }
        y[i] = ls;
    }
}

/* computes ||Ax-b||*/
double helper(const int n, const double *A, const double *b, const double *x, double *y){
    matrix_vector_mult(n, A, x, y);
    for (int i = 0; i < n; i++)
        y[i] -= b[i];

    double norm = 0.0;
    for (int i = 0; i < n; i++){
        norm += (y[i] * y[i]);
    }
    norm = sqrt(norm);
    return norm;
}
// implements the sequential jacobi method
void jacobi(const int n, double *A, double *b, double *x, int max_iter, double l2_termination)
{
    for (int i = 0; i < n; i++)
        x[i] = 0.0;
    double *D = (double*) (malloc(n * n * sizeof(double)));
    double *R = (double*) (malloc(n * n * sizeof(double)));
    double *invD = (double*) (malloc(n * n * sizeof(double)));
    memset(D, 0, n * n * sizeof(double));
    memset(R, 0, n * n * sizeof(double));
    memset(invD, 0, n * n * sizeof(double));
    
    for (int i = 0; i < n; i++)
        D[n * i + i] = A[n * i + i];
    for (int i = 0; i < n;i++){
        for (int j = 0; j < n;j++){
            int idx = n * i + j;
            R[idx] = A[idx] - D[idx];
        }
    }
    for (int i = 0; i < n; i++){
        if(D[n*i +i] == 0)
            return;
        invD[n * i + i] = 1 / D[n * i + i];
    }

    int itr = 0;
    double *y = (double *)(malloc(n * sizeof(double)));
    
    bool terminate = false;
    while(!terminate && (itr++) < max_iter){
        //y = Rx
        matrix_vector_mult(n, R, x, y);
        //if(itr==2)    printf("!!!!!!!!!!!!!!!!! %f\n",y[0]);
        //y = b - Rx
        for (int i = 0; i < n; i++)
            y[i] = b[i] - y[i];
        //x = (Dinv)y
        matrix_vector_mult(n, invD, y, x); 
        
        double norm = helper(n, A, b, x, y);
//         if(itr==3){
//             printf("SEQ NORM = %f %d\n",norm,itr);
//             for(int i=0;i<16;i++){
//                 printf("%f %d\n",x[i],i);
//             }
//         }
        if(norm <= l2_termination){
            terminate = true;
            //printf("SEQ NORM = %f %d\n",norm,itr);
        }
    }
    free(D); free(R); free(invD); free(y);
}
