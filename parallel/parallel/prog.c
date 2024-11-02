#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include <omp.h>
#include <string.h>
#include <time.h>
#include <ctype.h>
#include "cudaHeader.h"

enum constants
{
    MASTER = 0,
    MAX_SEQ1 = 3000,
    MAX_SEQ2 = 2000,
    RES_INFO = 3,
    MAT_SIZE = 26
};

int score_mat[MAT_SIZE][MAT_SIZE];

void create_default_matrix(int matrix[MAT_SIZE][MAT_SIZE])
{
    for (int i = 0; i < MAT_SIZE; i++)
    {
        for (int j = 0; j < MAT_SIZE; j++)
        {
            if (i == j)
            {
                matrix[i][j] = 1;
            }
            else
            {
                matrix[i][j] = 0;
            }
        }
    }
}

void print_matrix(int matrix[MAT_SIZE][MAT_SIZE])
{
    for (int i = 0; i < MAT_SIZE; i++)
    {
        for (int j = 0; j < MAT_SIZE; j++)
        {
            printf("%d ", matrix[i][j]);
        }
        printf("\n");
    }
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

void to_upper(char *seq1)
{
    for (int i = 0; i < strlen(seq1); i++)
    {
        seq1[i] = toupper(seq1[i]);
    }
}

int* omp_score_offset_mutant(int score_mat[MAT_SIZE][MAT_SIZE], char *seq1, char *seq2, int offset_max)
{
    int lenght_seq2 = strlen(seq2);
    int possible_mutants = lenght_seq2;
    int curr_score = 0;
    int index_seq1 = 0;

    int* res = (int*)malloc(sizeof(int) * RES_INFO);
    res[0] = INT_MIN;

#pragma omp parallel for shared(score_mat,seq1, seq2, lenght_seq2, possible_mutants) firstprivate(index_seq1, curr_score) 
    for (int offset = 0; offset < offset_max; offset++) //all offsets
    {
        for (int mutant = 0; mutant < possible_mutants; mutant++) //all mutants
        {
            index_seq1 = offset;
            char *new_seq = gen_mutant(seq2, mutant);
            for (int index_seq2 = 0; index_seq2 < lenght_seq2; index_seq2++, index_seq1++)
            {
                curr_score += score_mat[(seq1[index_seq1]) - 'A'][(new_seq[index_seq2]) - 'A'];
            }
            #pragma omp critical
            {
                if(curr_score > res[0])
                {
                    res[0] = curr_score;
                    res[1] = offset;
                    res[2] = mutant;
                }
            }
            curr_score = 0;
        }
    }
    return res;
}

void workerProcess(int score_mat[MAT_SIZE][MAT_SIZE], int num_of_seq, char* seq1)
{
    char* seq2 = (char*)malloc(sizeof(char) * MAX_SEQ2);
    int* res = (int*)malloc(sizeof(int) * RES_INFO);

    int* omp_res = (int*)calloc(sizeof(int) , RES_INFO);
    int* cuda_res = (int*)calloc(sizeof(int) , RES_INFO);
    
    int tag = 0;
    MPI_Status status;
    

    do
    {
        MPI_Recv(seq2,MAX_SEQ2,MPI_CHAR,MASTER,MPI_ANY_TAG,MPI_COMM_WORLD,&status); //RECEIVE WORK
        tag = status.MPI_TAG;
        
        if(tag != num_of_seq + 1) //NO SIGNAL TO DIE
        {
            //SPLIT WORK BETWEEN OMP AND CUDA
            int offsets = strlen(seq1) - strlen(seq2);
            int cuda_offset = offsets / 2 ; 
            int omp_offset = offsets / 2 + offsets % 2;

            //GET MAX SIZE FROM CUDA AND OMP
            omp_res = omp_score_offset_mutant(score_mat, seq1, seq2,omp_offset);
            cuda_res = cuda_score_offset_mutant(score_mat, seq1, seq2, omp_offset, offsets); //cuda off set start at number of omp offsets
            
            if(omp_res[0] < cuda_res[0])
                res = cuda_res;
            else
                res = omp_res;
            
            MPI_Send(&tag, 1, MPI_INT, MASTER, tag, MPI_COMM_WORLD); //SEND THE LOCATION
            MPI_Recv(&tag, 1, MPI_INT, MASTER, MPI_ANY_TAG, MPI_COMM_WORLD, &status); //GET APPROVED 
            MPI_Send(res, RES_INFO, MPI_INT, MASTER, tag, MPI_COMM_WORLD); //SEND THE DATA  
        }
        
    } while (tag != num_of_seq + 1);
    
    //FREE DATA
    free(res);
    free(seq2);
}


void masterProcess(int num_of_seq, char** all_seq, char* seq1, int num_proc)
{
    //SET THE NUMBERS
    int jobs_total = num_of_seq;
    int num_workers = num_proc - 1;
    int send_and_rec = 0, jobs_sent = 0;
    int worker_id;
    MPI_Status status;

    //GATHER ALL RESULTS OF OFFSET, MUTATION AND SCORE
    int** res = (int**)malloc(sizeof(int*) * num_of_seq);

    for (int i = 0; i < num_of_seq; i++)
    {
        res[i] = (int*)malloc(sizeof(int) * RES_INFO);
    }
    

    //START WORKERS
    for(worker_id = 1; worker_id < num_proc; worker_id++)
    {
        MPI_Send(all_seq[jobs_sent], MAX_SEQ2, MPI_CHAR, worker_id, jobs_sent, MPI_COMM_WORLD);
        jobs_sent ++;
        send_and_rec ++;
    }
    
    //RECIVE AND SEND MORE WORK
    while(send_and_rec != 0)
    {
        //GET THE DATA AND SEND MORE
        int location = 0;
        MPI_Recv(&location, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status); //GET THE LOCATION
        MPI_Send(&location, 1, MPI_INT, status.MPI_SOURCE, location, MPI_COMM_WORLD);
        MPI_Recv(res[location], RES_INFO, MPI_INT, status.MPI_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status); //GET THE DATA
        send_and_rec --;
        
        if((jobs_sent < jobs_total))
        {
            MPI_Send(all_seq[jobs_sent], MAX_SEQ2, MPI_CHAR, status.MPI_SOURCE, jobs_sent, MPI_COMM_WORLD);
            jobs_sent ++;
            send_and_rec ++;
        }
    }
    
    //SEND THE SIGNAL TO KILL PROCESSES
    int stop = num_of_seq + 1;
    for (worker_id = 1; worker_id < num_proc; worker_id++)
    {
        MPI_Send(&stop, 1 ,MPI_INT, worker_id, stop, MPI_COMM_WORLD);
    }
    
    //PRINT THE RESULT
    
    printf("-----------------------------------Parallel Program-----------------------------------------------\n");
    printf("\n");
    printf("Seq1 is %s\n", seq1);
    printf("\nHere are all the best matches for the strings:\n");
    
    for (int i = 0; i < num_of_seq; i++)
    {
        printf("seq2 = %s\n", all_seq[i]);
        if(res[i][2] == 0) //if mutation = 0 then change it to max length
            res[i][2] = strlen(all_seq[i]);
        printf("offset (n) = %d , mutation (k) = %d , score = %d \n", res[i][1], res[i][2], res[i][0]);
        printf("\n");
    }
    
    //FREE DATA
    for (int i = 0; i < num_of_seq; i++)
    {
        free(res[i]);
    }
    free(res);
}





int main(int argc, char *argv[])
{
    double start = time(NULL);
    //GENERAL DATA NEEDED
    int my_rank, num_procs, number_of_sequences;
    char* seq1 = (char*)malloc(sizeof(char) * MAX_SEQ1);
    char** all_seq;

    
    //MPI INIT
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs); 


    if(my_rank == MASTER)
    {
        //GET ALL THE INPUT

        //GET SEQ1
        scanf("%s", seq1);

        //GET NUMBER OF SEQ'S
        scanf("%d", &number_of_sequences);

        //GET ALL OTHER SEQ'S
            //MALLOC **
        all_seq = (char**)malloc(sizeof(char*) * number_of_sequences);
            //MALOC *
        for (int i = 0; i < number_of_sequences; i++)
        {
            all_seq[i] = (char*)malloc(sizeof(char) * MAX_SEQ2);
        }
            //GET FROM FILE
        for (int i = 0; i < number_of_sequences; i++)
        {
            scanf("%s",all_seq[i]);
        }

        // CHANGE SEQs TO UPPER CASE
    to_upper(seq1);

    for (int i = 0; i < number_of_sequences; i++)
    {
        to_upper(all_seq[i]);
    }

    // GET THE SCORE TABLE
    if (argc != 2)
    {
        printf("Using the default score table\n");
        // CREATE DEFAULT MATRIX
        create_default_matrix(score_mat);
        // print_matrix(score_mat);
    }
    else
    {
        printf("Using custom score table\n");
        // OPEN FILE
        FILE *fp = fopen(argv[1], "r");
        if (fp == NULL)
        {
            printf("Error opening file\n");
            exit(1);
        }
        // READ FROM FILE
        for (int i = 0; i < MAT_SIZE; i++)
        {
            for (int j = 0; j < MAT_SIZE; j++)
            {
                fscanf(fp, "%d", &score_mat[i][j]);
            }
        }
        // CLOSE FILE
        fclose(fp);
    }
        
        //MAKE SURE DATA IS BIG ENOUGH
         if(number_of_sequences < num_procs - 1)
        {
            printf("The program needs at least the %d inputs to work\n", num_procs - 1);
            MPI_Finalize();
            return 0;
        }
    }
    
   

    //BROADCAST SEQ1 AND SCORE MATRIX TO ALL PROCESSES
    MPI_Bcast(seq1, MAX_SEQ1, MPI_CHAR, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(score_mat, MAT_SIZE*MAT_SIZE, MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(&number_of_sequences, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
    
    //START THE LOGIC --------------------> MASTER WORKER DYNAMIC ARCHITECTURE 
    if(my_rank == MASTER)
    {
        masterProcess(number_of_sequences, all_seq, seq1, num_procs);
    }
    else
    {
        workerProcess(score_mat, number_of_sequences, seq1);
    }
    
    //FREE ALL DATA
    
    if(my_rank == MASTER)
    {
        for (int i = 0; i < number_of_sequences; i++)
        {
            free(all_seq[i]);
        }
        free(all_seq);
    }

    free(seq1);
    
    if(my_rank == MASTER)
        printf("the time it took is %lf seconds\n", time(NULL) - start);
    
    MPI_Finalize();
    return 0;
}
