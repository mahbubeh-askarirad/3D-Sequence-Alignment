#ifndef DYNAMIC_PROGRAMMING_H
#define DYNAMIC_PROGRAMMING_H
#include "Alignment_params.h"
#include <stdio.h>

void dynamic_programming(int* s1, int* s2, int* s3, int seq1_length, 
int seq2_length, int seq3_length, int*** AlignedTriplet, AlignmentParams params, int* triplet_length);

#endif