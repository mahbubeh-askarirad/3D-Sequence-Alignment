/*
 * Triple_Alignment.c
 * Multiple sequence alignment for 3 sequences with affine gap penalties, 
 * using three dimentional Dynamic Programing 
 * Author: Mahbubeh Askarirad
 * Date: 2024-06-28
 * License: MIT
 *
 * Usage: Compile with all required modules (see README.md)
 *        Reads input sequences from a FASTA-like file.
 *        Outputs alignment and scores to alignment_output.txt.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Readsequences.h"
#include "write_sequence.h"
#include "mmMatrix.h"
#include "Alignment_params.h"


#define MAX_SEQ_LENGTH 1000
#define NUM_SEQS 3

void Readsequences(FILE* output, const char* filename, int*** Seqs, int** lengths);
void write_sequence(FILE* output, int l, int* seq );
void Alignment_scoring_matrix(FILE* output, int* seq1, int* seq2, int* rev_seq1, int* rev_seq2, 
    int seq1_length, int seq2_length, int*** Aligned_pair, int*** score_matrix, AlignmentParams params, int* score, int* alignment_2D_l);
void dynamic_programming(int* s1, int* s2, int* s3, int seq1_length, int seq2_length, int seq3_length, int*** AlignedTriplet, AlignmentParams params, int* triplet_length);
int sum_of_pair(int* s1, int* s2, int* s3, int l, AlignmentParams params);


int Scoring_3d(int** sequences, int length, AlignmentParams params){
    int score = 0;
    int s = 0;
    int* s1 = sequences[0];
    int* s2 = sequences[1];
    int* s3 = sequences[2];

    if(s1[0] != 0 && s2[0] != 0 && s3[0] != 0) // M
        score = mmMatrix(params.mmflag, s1[0], s2[0]) + 
            mmMatrix(params.mmflag, s1[0], s3[0])+
            mmMatrix(params.mmflag, s2[0], s3[0]);

    else if(s1[0] != 0 && s2[0] != 0 && s3[0] == 0) //XY
        score = mmMatrix(params.mmflag, s1[0], s2[0]) - 2 * params.GO;

    else if(s1[0] != 0 && s2[0] == 0 && s3[0] != 0) //XZ
        score = mmMatrix(params.mmflag, s1[0], s3[0]) - 2 * params.GO;

    else if(s1[0] == 0 && s2[0] != 0 && s3[0] != 0) //YZ
        score = mmMatrix(params.mmflag, s2[0], s3[0]) - 2 * params.GO;

    else
        score = - 2 * params.GO;

    for(int i = 1; i < length; i++){ 
        
        s = 0;
        
        if(s1[i] != 0 && s2[i] != 0 && s3[i] != 0){ // XYZ
            
            s = mmMatrix(params.mmflag, s1[i], s2[i]) + 
                mmMatrix(params.mmflag, s1[i], s3[i])+
                mmMatrix(params.mmflag, s2[i], s3[i]);
        }
        if(s1[i] != 0 && s2[i] != 0 && s3[i] == 0){ //XY
            if(s3[i-1] == 0){
                s = mmMatrix(params.mmflag, s1[i], s2[i]) - 2 * params.GE;
            }
            else{
                s = mmMatrix(params.mmflag, s1[i], s2[i]) - 2 * params.GO;
            }
        }
        else if(s1[i] != 0 && s2[i] == 0 && s3[i] != 0){//XZ
            if(s2[i-1] == 0){
                s = mmMatrix(params.mmflag, s1[i], s3[i]) - 2 * params.GE;
            }
            else{
                s = mmMatrix(params.mmflag, s1[i], s3[i]) - 2 * params.GO;
            }
        }
        else if(s1[i] == 0 && s2[i] != 0 && s3[i] != 0){//YZ
            if(s1[i-1] == 0){
                s = mmMatrix(params.mmflag, s2[i], s3[i]) - 2 * params.GE;
            }
            else{
                s = mmMatrix(params.mmflag, s2[i], s3[i]) - 2 * params.GO;
            }
        }
        else if(s1[i] != 0 && s2[i] == 0 && s3[i] == 0){//X
            if(s2[i-1] == 0 && s3[i-1] == 0){
                s = - 2 * params.GE;
            }
            else if(s2[i-1] != 0 && s3[i-1] != 0){
                s = - 2 * params.GO;
            }
            else{
                s = - params.GE - params.GO;
            }
        }
        else if(s1[i] == 0 && s2[i] != 0 && s3[i] == 0){//Y
            if(s1[i-1] == 0 && s3[i-1] == 0){
                s = - 2 * params.GE;
            }
            else if(s1[i-1] != 0 && s3[i-1] != 0){
                s = - 2 * params.GO;
            }
            else{ 
                s = - params.GE - params.GO;
            }
        }
        else if(s1[i] == 0 && s2[i] == 0 && s3[i] != 0){//Z
            if(s1[i-1] == 0 && s2[i-1] == 0){
                s = - 2 * params.GE;
            }
            else if(s1[i-1] != 0 && s2[i-1] != 0){
                s = - 2 * params.GO;
            }
            else{
                s = - params.GE - params.GO;
            }
        }
        score = score +  s;
    }
    return score;
}

int sum_of_pair(int* seq_1, int* seq_2, int* seq_3, int l, AlignmentParams params){

    int score = 0;
    int sum = 0;
    int gflag1 = 1;
    int gflag2 = 1;

    int* st1 = NULL;
    int* st2 = NULL;

    int* s1;
    int* s2;

    s1 = malloc((sizeof(int)) * l);
    s2 = malloc((sizeof(int)) * l);
    
    for(int pair = 0; pair < 3; pair++){

        //printf("pair : %d\n\n", pair);
        
        switch (pair)
        {
        case 0:
            st1 = seq_1;
            st2 = seq_2;
            break;
        case 1:
            st1 = seq_1;
            st2 = seq_3;
            break;
        case 2:
            st1 = seq_2;
            st2 = seq_3;
            break;
        
        default:
            break;
        }

        int ctr = 0;
        for(int p = l-1; p >= 0; p--){
            if(st1[p] != 0 || st2[p] != 0){
                s1[ctr] = st1[p];
                s2[ctr] = st2[p];
                ctr++;
            }
        }
        int ll = ctr;
        score = 0;

        for(int i=0; i<ll; i++){

            if(s1[i] != 0 && s2[i] != 0){ 
                score = score + mmMatrix( params.mmflag, s1[i], s2[i]);
                gflag1 = 1;
                gflag2 = 1;
            }else if(s1[i] == 0 && s2[i] != 0){
                gflag2 = 1;
                if(gflag1 == 1){
                    gflag1 = 0;
                    score = score - params.GO;
                }
                else
                    score = score - params.GE;
            } else if(s1[i] != 0 && s2[i] == 0){
                gflag1 = 1;
                if(gflag2 == 1){
                    gflag2 = 0;
                    score = score - params.GO;
                }
                else
                    score = score - params.GE;
            }
        }

        sum = sum + score;  
    }
    return sum;
}

int read_positive_int(const char* prompt){

    char buffer[100];
    int value = -1;
    while (value < 0){
        printf("%s", prompt);
        if(fgets(buffer, sizeof(buffer), stdin)){
            if(sscanf(buffer, "%d", &value) == 1 && value < 0) {
                printf("Invalid input. Please enter a non-negative integer.\n");
                value = -1;
            }
        }
    }
    return value;
}

int select_matrix(){
    char buffer[100];
    int choice = -1;
    while (choice != 0 && choice != 1 && choice != 2){
        printf("select scoring matrix:\n");
        printf(" 0: BLOSUM62\n");
        printf(" 1: BLOSUM30\n");
        printf(" 2: BLOSUM90\n");
        printf("Enter your choice: ");

        if(fgets(buffer, sizeof(buffer), stdin)){
            if(sscanf(buffer, "%d", &choice) != 1 || 
              (choice != 0 && choice != 1 && choice != 2)){
                printf("Invalid choice. Please enter: 0, 1 or 2.\n");
                choice = -1;
              }
        }
    }
    return choice;

}

AlignmentParams get_alignment_params(){
    AlignmentParams params;
    params.GO = read_positive_int("Enter gap opening penalty (GO): ");
    params.GE = read_positive_int("Enter gap extension penalty (GE): ");
    params.mmflag = select_matrix();
    return params;
}

int main(){

    int** Seqs = NULL;
    int* seq_lengths = NULL;
    int** AlignedTriplet = NULL;
    int triplet_length = 0;

    FILE* output = fopen("alignment_output.txt", "w");
    if (output == NULL) {
        fprintf(output, "Failed to open output file.\n");
        return 1;
    }

    AlignmentParams params = get_alignment_params();
    printf("\n--- Alignemnt Settings ---\n");
    printf("Gap Opening Penalty: %d\n", params.GO);
    printf("Gap Extension Penalty: %d\n", params.GE);
    printf("Scoring Matrix: %s", 
        params.mmflag == 0 ? "BLOSUM62\n":
        params.mmflag == 1 ? "BLOSUM30\n":  "BLOSUM90\n");
    printf("\n");
    
    Readsequences(output, "input", &Seqs, &seq_lengths);  // Reading sequences from FASTA file

    for (int i = 0; i < NUM_SEQS; i++) {  // seqeunce length
        fprintf(output, "Sequence %d length: %d\n\n", i + 1, seq_lengths[i]);
    }
    int seq1_length = seq_lengths[0];
    int seq2_length = seq_lengths[1];
    int seq3_length = seq_lengths[2];

    dynamic_programming(Seqs[0], Seqs[1], Seqs[2], seq1_length, seq2_length, seq3_length, 
                        &AlignedTriplet, params, &triplet_length);

    // reverse sequences
    int a, b ,c;
    for(int i =0; i < triplet_length / 2; i++){
        a = AlignedTriplet[0][i];
        b = AlignedTriplet[1][i];
        c = AlignedTriplet[2][i];

        AlignedTriplet[0][i] = AlignedTriplet[0][triplet_length - (i + 1)];
        AlignedTriplet[1][i] = AlignedTriplet[1][triplet_length - (i + 1)];
        AlignedTriplet[2][i] = AlignedTriplet[2][triplet_length - (i + 1)];

        AlignedTriplet[0][triplet_length - (i + 1)] = a;
        AlignedTriplet[1][triplet_length - (i + 1)] = b;
        AlignedTriplet[2][triplet_length - (i + 1)] = c;
    }
    int triplet_3d_score = Scoring_3d(AlignedTriplet, triplet_length, params);
    int triplet_alignment_score = sum_of_pair(AlignedTriplet[0], AlignedTriplet[1], AlignedTriplet[2], triplet_length, params);

    fprintf(output, "Alignment length:  %d\n", triplet_length);
    fprintf(output, "\nAlignment SPS : %d\n", triplet_alignment_score);
    fprintf(output,"\nAlignment TPS %d\n", triplet_3d_score );
    
    fprintf(output, "\n\n\n>sequence_1\n");
    write_sequence(output, triplet_length, AlignedTriplet[0]);
    fprintf(output, ">sequence_2\n");
    write_sequence(output, triplet_length, AlignedTriplet[1]);
    fprintf(output, ">sequence_3\n");
    write_sequence(output, triplet_length, AlignedTriplet[2]);

    fclose(output);

    // ------- FREE ARRAYS

    for (int i = 0; i < NUM_SEQS; i++) {
    free(AlignedTriplet[i]);
    }
    free(AlignedTriplet);

    // Free Seqs (the input sequences)
    for (int i = 0; i < NUM_SEQS; i++) {
        free(Seqs[i]);
    }
    free(Seqs);

    // Free sequence lengths
    free(seq_lengths);
    return 0;
}


