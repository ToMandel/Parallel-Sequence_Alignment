#pragma once

#define NUM_THREADS_PER_BLOCK 256

const int CUDA_MAT_SIZE = 26;

int* cuda_score_offset_mutant(int score_mat[CUDA_MAT_SIZE][CUDA_MAT_SIZE], char* seq1, char* seq2, int offset_start, int offset_finish);
