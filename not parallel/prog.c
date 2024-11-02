#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <ctype.h>
#include <time.h>

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
        if (!new_seq)
        {
            printf("Error allocating memory\n");
            exit(1);
        }
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

int *score_offset_mutant(int score_mat[MAT_SIZE][MAT_SIZE], char *seq1, char *seq2, int offset_max)
{
    int lenght_seq2 = strlen(seq2);
    int possible_mutants = lenght_seq2;
    int curr_score = 0;
    int index_seq1 = 0;

    int *res = (int *)malloc(sizeof(int) * RES_INFO);
    if (!res)
    {
        printf("Error allocating memory\n");
        exit(1);
    }
    res[0] = INT_MIN;

    for (int offset = 0; offset <= offset_max; offset++) // all offsets
    {
        for (int mutant = 0; mutant <= possible_mutants; mutant++) // all mutants
        {
            index_seq1 = offset;
            char *new_seq = gen_mutant(seq2, mutant);
            for (int index_seq2 = 0; index_seq2 < lenght_seq2; index_seq2++, index_seq1++)
            {
                curr_score += score_mat[(seq1[index_seq1]) - 'A'][(new_seq[index_seq2]) - 'A'];
            }
            if (curr_score > res[0])
            {
                res[0] = curr_score;
                res[1] = offset;
                res[2] = mutant;
            }
            curr_score = 0;
            free(new_seq);
        }
    }
    return res;
}

void to_upper(char *seq)
{
    for (int i = 0; i < strlen(seq); i++)
    {
        seq[i] = toupper(seq[i]);
    }
}

int main(int argc, char *argv[])
{
    double start = time(NULL);

    // GENERAL DATA NEEDED
    int number_of_sequences;
    char **all_seq;

    // GET ALL THE INPUT
    // GET SEQ1
    char seq1[MAX_SEQ1];
    scanf("%s", seq1);

    // GET NUMBER OF SEQ'S
    scanf("%d", &number_of_sequences);

    // GET ALL OTHER SEQ'S
    // MALLOC **
    all_seq = (char **)malloc(sizeof(char *) * number_of_sequences);
        if (!all_seq)
        {
            printf("Error allocating memory\n");
            exit(1);
        }
    // MALOC *
    for (int i = 0; i < number_of_sequences; i++)
    {
        all_seq[i] = (char *)malloc(sizeof(char) * MAX_SEQ2);
                if (!all_seq[i])
        {
            printf("Error allocating memory\n");
            exit(1);
        }
    }
    // GET FROM FILE
    for (int i = 0; i < number_of_sequences; i++)
    {
        scanf("%s", all_seq[i]);
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

    // RESULTS
    int **res = (int **)malloc(sizeof(int *) * number_of_sequences);
            if (!res)
        {
            printf("Error allocating memory\n");
            exit(1);
        }
    
    // IN CASE THERE IS MORE THAN ONE SEQUENCE
    for (int i = 0; i < number_of_sequences; i++)
    {
        res[i] = (int *)malloc(sizeof(int) * RES_INFO);
        if (!res[i])
        {
            printf("Error allocating memory\n");
            exit(1);
        }
    }

    for (int i = 0; i < number_of_sequences; i++)
    {
        int offset_max = strlen(seq1) - strlen(all_seq[i]);
        res[i] = score_offset_mutant(score_mat, seq1, all_seq[i], offset_max);
    }

    // PRINT THE RES
    printf("-----------------------------------Non Parallel Program-----------------------------------------------\n");
    printf("\n");
    printf("Seq1 is %s\n", seq1);
    printf("\nHere are all the best matches for the strings:\n");

    for (int i = 0; i < number_of_sequences; i++)
    {
        printf("seq2 = %s\n", all_seq[i]);
        printf("highest alignment score = %d , offset = %d , k = %d \n", res[i][0], res[i][1], res[i][2]);
        printf("\n");
    }

    // FREE ALL DATA
    printf("Freeing all_seq\n");
    free(all_seq);

    printf("Freeing res\n");
    free(res);

    printf("elapsed time is %lf seconds\n", time(NULL) - start);
    return 0;
}