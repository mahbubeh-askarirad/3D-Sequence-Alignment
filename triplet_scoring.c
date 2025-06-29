#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Readsequences_with_gaps.h"
#include "write_sequence.h"
#include "mmMatrix.h"
#include "PairDistance.h"
#include "lambda_e_estimation.h"
#include "Alignment_scoring_matrix.h"
#include "Alignment_params.h"


#define MAX_SEQ_LENGTH 1000


int Scoring_3d(int** sequences, int length, AlignmentParams params){
    int score = 0;
    int s = 0;
    int* s1 = sequences[0];
    int* s2 = sequences[1];
    int* s3 = sequences[2];
    //printf("%d. %d. %d\n\n", s1[length-1], s2[length-1], s3[length-1]);

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
            else{ printf("HELLO\n");
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
        //printf("%d    %d %d %d\n", s, s1[i],s2[i],s3[i]);
        score = score +  s;
        printf("%d  %d    %d %d %d\n", s, score, s1[i],s2[i],s3[i]);
    }
    return score;
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

AlignmentParams get_alignemnt_params(){
    AlignmentParams params;
    params.GO = read_positive_int("Enter gap opening penalty (GO): ");
    params.GE = read_positive_int("Enter gap extension penalty (GE): ");
    params.mmflag = select_matrix();
    return params;
}

int main (){

    printf("ERROR\n");

    AlignmentParams params = get_alignemnt_params();

    int** Seqs = NULL;
    int** rev_Seqs = NULL;
    int* seq_lengths = NULL;
    int num_seqs = 0;
    int** Aligned_pair = NULL;

    Readsequences_with_gaps("MAFFT1", &Seqs, &seq_lengths, &num_seqs);  // Reading sequences from FASTA file

    int score = Scoring_3d(Seqs, seq_lengths[0], params);

    printf("\ntriplet_score: %d\n", score);

}