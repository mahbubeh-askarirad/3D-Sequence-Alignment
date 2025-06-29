#ifndef ALIGNMENT_SCORING_MATRIX_H
#define ALIGNMENT_SCORING_MATRIX_H
#include <stdio.h>
#include "Alignment_params.h"

void Alignment_scoring_matrix(FILE* output, int* seq1, int* seq2, int* rev_seq1, int* rev_seq2, 
    int seq1_length, int seq2_length, int*** Aligned_pair, int*** score_matrix, AlignmentParams params, int* score, int* alignment_2D_l);

#endif