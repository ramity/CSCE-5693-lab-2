/******************************************************************************
 *cr
 *cr            (C) Copyright 2010 The Board of Trustees of the
 *cr                        University of Illinois
 *cr                         All Rights Reserved
 *cr
 ******************************************************************************/

#include <stdio.h>

__global__ void mysgemm(int m, int n, int k, const float *A, const float *B, float* C)
{
    /********************************************************************
     *
     * Compute C = A x B
     *   where A is a (m x k) matrix
     *   where B is a (k x n) matrix
     *   where C is a (m x n) matrix
     *
     ********************************************************************/

    // INSERT KERNEL CODE HERE
    int row = blockIdx.y * blockDim.y + threadIdx.y;
    int col = blockIdx.x * blockDim.x + threadIdx.x;
    int size = m * n;

    if(row < m && col < n)
    {
        int sum = 0;

        for(int i = 0;i < size;i++)
        {
            sum += A[row * k + i] * B[i * n + col];
        }

        C[row * n + col] = sum;
    }
}

void basicSgemm(char transa, char transb, int m, int n, int k, float alpha, const float *A, int lda, const float *B, int ldb, float beta, float *C, int ldc)
{
    if ((transa != 'N') && (transa != 'n'))
    {
        printf("unsupported value of 'transa'\n");
        return;
    }

    if ((transb != 'N') && (transb != 'n'))
    {
        printf("unsupported value of 'transb'\n");
        return;
    }

    if ((alpha - 1.0f > 1e-10) || (alpha - 1.0f < -1e-10))
    {
        printf("unsupported value of alpha\n");
        return;
    }

    if ((beta - 0.0f > 1e-10) || (beta - 0.0f < -1e-10))
    {
        printf("unsupported value of beta\n");
        return;
    }

    // Initialize thread block and kernel grid dimensions ---------------------

    const unsigned int BLOCK_SIZE = 16; // Use 16x16 thread blocks

    //INSERT CODE HERE - UNSURE IF CORRECT; REVIEW IF VERIFY ERRORS OCCUR
    int C_size = (m * n);
    dim3 blockDim(BLOCK_SIZE, BLOCK_SIZE);
    int numberOfBlocks = ceil(C_size / BLOCK_SIZE);
    dim3 gridDim(numberOfBlocks, numberOfBlocks);

    // Invoke CUDA kernel -----------------------------------------------------

    //INSERT CODE HERE
    mysgemm<<<gridDim,blockDim>>>(m, n, k, A, B, C);

}
