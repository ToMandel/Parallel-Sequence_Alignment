#include <cuda_runtime.h>
// #include <helper_cuda.h>
#include <stdio.h>
#include <iostream>
#include <string.h>
#include <limits.h>
#include "cudaHeader.h"

enum constants
{
    RES_INFO = 3,
    MAT_SIZE = 26,
    MAX_SEQ2 = 2000,
};

__device__ int score_mat[MAT_SIZE][MAT_SIZE];
__device__ __constant__ int res_components = 3;

__device__ int cuda_strlen(char* seq)
{
    int counter = 0;
    while(*seq != '\0')
    {
        counter++;
        seq++;
    }

    return counter;
}
char *gen_mutant(char *seq2, int mutant)
{
    char *new_seq = (char *)malloc(sizeof(char) * MAX_SEQ2);
    strcpy(new_seq, seq2);
    for (int i = mutant; i < strlen(seq2); i++)
    {
        if (new_seq[i] >= 'Z')
        {
            new_seq[i] = 'A';
        }
        else
        {
            new_seq[i] = new_seq[i] + 1;
        }
    }
    return new_seq;
}

__global__ void score_offset_mutant(int score_mat[MAT_SIZE][MAT_SIZE],char* seq1, char* seq2 , int* lenght_seq2 , int* start, int* finish, int* res)
{
    int possible_mutants = *lenght_seq2;
    int curr_score = 0;
    int index = 0;

    int thread_id = threadIdx.x + blockIdx.x * blockDim.x;
    int offset =  thread_id + *start;
    int thread_index_for_res_array = thread_id * res_components; 

    if(offset <= *finish)
    {
        for (int mutant = 0; mutant < possible_mutants; mutant++) //all mutants
        {
            index = offset;
            char *new_seq = gen_mutant(seq2, mutant);
            for (int index_seq2 = 0; index_seq2 < *lenght_seq2; index_seq2++, index++)
            {
                curr_score += score_mat[(seq1[index]) - 'A'][(new_seq[index_seq2]) - 'A'];
            }
            
            if(curr_score > res[thread_index_for_res_array])
            {
                res[thread_index_for_res_array] = curr_score;
                res[thread_index_for_res_array + 1] = offset;
                res[thread_index_for_res_array + 2] = mutant;
            }
            curr_score = 0;
        }
    } 
}

int* cuda_score_offset_mutant(int score_mat[MAT_SIZE][MAT_SIZE],char* seq1, char* seq2 ,int offset_start, int offset_finish)
{
    //CREATE THE ARRAY THAT WILL RETURN THE BEST RESULT
    int* result = (int*)calloc(sizeof(int), RES_INFO);

    //DATA NEEDED
    int offset_size = offset_finish - offset_start + 1; //--> num of cuda threads
    int num_of_blocks = (offset_size / NUM_THREADS_PER_BLOCK);
    if (offset_size % NUM_THREADS_PER_BLOCK != 0)
        num_of_blocks ++;
    
    int lenght_of_res = (RES_INFO) * offset_size; //for each offset the gpu finds the best mutant

    //ALLOCATE DATA TO CUDA MEMORY
    char* cuda_seq1, *cuda_seq2;
    int cuda_matrix[MAT_SIZE][MAT_SIZE];
    int* cuda_offset_start, *cuda_offset_finish ,*cuda_seq2_lenght;
    int* cuda_res, *res = (int*)calloc(sizeof(int) , lenght_of_res);
    int seq2_lenght = (strlen(seq2));

    for (int i = 0; i < lenght_of_res; i++) // get the array ready
    {
        res[i] = INT_MIN;
    }

        //sizes to allocate
    int size_for_cuda_seq1 = sizeof(char) * (strlen(seq1));
    int size_for_cuda_seq2 = sizeof(char) * seq2_lenght;
    int size_for_cuda_matrix = sizeof(int) * (MAT_SIZE*MAT_SIZE);
    int size_for_cuda_int = sizeof(int);
    int size_for_cuda_res = sizeof(int) * lenght_of_res;

        //allocate
    cudaMalloc((void**)&cuda_seq1, size_for_cuda_seq1);
    cudaMalloc((void**)&cuda_seq2, size_for_cuda_seq2);
    cudaMalloc((void**)&cuda_offset_start, size_for_cuda_int);
    cudaMalloc((void**)&cuda_offset_finish, size_for_cuda_int);
    cudaMalloc((void**)&cuda_res, size_for_cuda_res);
    cudaMalloc((void**)&cuda_seq2_lenght, size_for_cuda_int);



    //COPY INPUT INTO DEVICE
    cudaMemcpy(cuda_seq1, seq1, size_for_cuda_seq1, cudaMemcpyHostToDevice);
    cudaMemcpy(cuda_seq2, seq2, size_for_cuda_seq2, cudaMemcpyHostToDevice);
    cudaMemcpy(cuda_matrix, score_mat, size_for_cuda_matrix, cudaMemcpyHostToDevice);
    cudaMemcpy(cuda_offset_start, &offset_start, size_for_cuda_int, cudaMemcpyHostToDevice);
    cudaMemcpy(cuda_offset_finish, &offset_finish, size_for_cuda_int, cudaMemcpyHostToDevice);
    cudaMemcpy(cuda_seq2_lenght, &seq2_lenght, size_for_cuda_int, cudaMemcpyHostToDevice);
    cudaMemcpy(cuda_res, res, size_for_cuda_res, cudaMemcpyHostToDevice);
    
    
    //LUNCH KERNEL
    score_offset_mutant<<<num_of_blocks, NUM_THREADS_PER_BLOCK>>>(cuda_seq1, cuda_seq2, cuda_seq2_lenght, **cuda_matrix, cuda_offset_start, cuda_offset_finish ,cuda_res);

    //COPY RESULT BACK TO HOST
    cudaMemcpy(res, cuda_res, size_for_cuda_res, cudaMemcpyDeviceToHost);
    
    //GET THE BIGGEST SCORE
    result[0] = INT_MIN; 

    for (int i = 0; i < lenght_of_res; i += 3)
    {
        if(result[0] < res[i])
        {
            result[0] = res[i];
            result[1] = res[i+1];
            result[2] = res[i+2];
        }
    }
    
    //FREE
    cudaFree(cuda_seq1);
    cudaFree(cuda_seq2);
    cudaFree(cuda_matrix);
    cudaFree(cuda_offset_start);
    cudaFree(cuda_offset_finish);
    cudaFree(cuda_res);
    cudaFree(cuda_seq2_lenght);

    return result;
}
