#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Alignment_scoring_matrix.h"
#include "write_sequence.h"
#include "mmMatrix.h"
#include "Alignment_params.h"

int max_two_number(int num1, int num2){
    if (num1 >= num2)
        return num1;
    else
        return num2;
}
int B0_func(int** B0, int** B1, int** B2, int*** Aligned_pair, int *s1, int *s2, int C, int R, int ct1, int ct2, int* final_length);
int B1_func(int** B0, int** B1, int** B2, int*** Aligned_pair, int *s1, int *s2, int C, int R, int ct1, int ct2, int* final_length);
int B2_func(int** B0, int** B1, int** B2,int*** Aligned_pair, int *s1, int *s2, int C, int R, int ct1, int ct2, int* final_length);

int B0_func(int** B0, int** B1, int** B2,int*** Aligned_pair, int *s1, int *s2, int C, int R, int ct1, int ct2, int* final_length){
    
    if(R == 0 && C == 0 || R < 0 || C < 0)
        return 0;

    (*Aligned_pair)[0][ct1] = s1[C-1];
    (*Aligned_pair)[1][ct2] = s2[R-1];
    ct1++;
    ct2++;
    (*final_length)++;
    
    if (B0[R][C] == 0){ 
        B0_func(B0, B1, B2, Aligned_pair, s1, s2, C-1, R-1, ct1, ct2, final_length);
    }
    else if(B0[R][C] == 1 || B0[R][C] == 3){
        B1_func(B0, B1, B2, Aligned_pair, s1, s2, C-1, R-1, ct1, ct2, final_length);
    }
    else if(B0[R][C] == 2){
        B2_func(B0, B1, B2, Aligned_pair, s1, s2, C-1, R-1, ct1, ct2, final_length);
    }
    return 0;

}

int B1_func(int** B0, int** B1, int** B2, int*** Aligned_pair, int *s1, int *s2, int C, int R, int ct1, int ct2, int* final_length){

    //printf("B1 %d   %d  %d\n", B1[R][C], R, C);
    if(R == 0 && C == 0 || R < 0 || C < 0)
        return 0;

    (*Aligned_pair)[0][ct1] = 0;
    (*Aligned_pair)[1][ct2] = s2[R-1];
    ct1++;
    ct2++;
    (*final_length)++;

    if (B1[R][C] == 0){
        
        B0_func(B0, B1, B2, Aligned_pair, s1, s2, C, R-1, ct1, ct2,final_length); 
    }
    else if(B1[R][C] == 1){
        B1_func(B0, B1, B2, Aligned_pair, s1, s2, C, R-1, ct1, ct2, final_length);
    }
    else    
        B2_func(B0, B1, B2, Aligned_pair, s1, s2, C, R-1, ct1, ct2, final_length);
    return 0;
}

int B2_func(int** B0, int** B1, int** B2, int*** Aligned_pair, int *s1, int *s2, int C, int R, int ct1, int ct2, int* final_length){

    //printf("B2 %d   %d  %d\n", B2[R][C], R, C);
    if(R == 0 && C == 0 || R < 0 || C < 0)
        return 0;

    (*Aligned_pair)[0][ct1] = s1[C-1];
    (*Aligned_pair)[1][ct2] = 0;
    ct1++;
    ct2++;
    (*final_length)++;

    if (B2[R][C] == 0){
        
        B0_func(B0, B1, B2, Aligned_pair, s1, s2, C-1, R, ct1, ct2, final_length); 
    }
    else if(B2[R][C] == 1){
        B1_func(B0, B1, B2, Aligned_pair, s1, s2, C-1, R, ct1, ct2, final_length);
    }
    else    
        B2_func(B0, B1, B2, Aligned_pair, s1, s2, C-1, R, ct1, ct2, final_length);

    return 0;
}

void Alignment_2D(int*** Aligned_pair, int*** M_score, int*** X_score, int*** Y_score, 
    int* s1, int* s2, int l1, int l2, int* final_length, AlignmentParams params){
        
    *Aligned_pair = malloc(2 * sizeof(int*));  // Two sequences in the pair
    (*Aligned_pair)[0] = malloc((l1 + l2) * sizeof(int));  // Worst case: all gaps
    (*Aligned_pair)[1] = malloc((l1 + l2) * sizeof(int));
   
    int ct1 = 0;
    int ct2 = 0;

    int** MM;
    int** X2;
    int** Y2;
    int** B0;
    int** B1;
    int** B2;

    MM = malloc(sizeof(int*) * (l2+1));
    X2 = malloc(sizeof(int*) * (l2+1));
    Y2 = malloc(sizeof(int*) * (l2+1));
    B0 = malloc(sizeof(int*) * (l2+1));
    B1 = malloc(sizeof(int*) * (l2+1));
    B2 = malloc(sizeof(int*) * (l2+1));
    *M_score = malloc(sizeof(int*) * (l2+1));

    for (int i = 0; i < l2+1; i++){
        MM[i] = malloc(sizeof(int) * (l1 + 1));
        X2[i] = malloc(sizeof(int) * (l1 + 1));
        Y2[i] = malloc(sizeof(int) * (l1 + 1));
        B0[i] = malloc(sizeof(int) * (l1 + 1));
        B1[i] = malloc(sizeof(int) * (l1 + 1));
        B2[i] = malloc(sizeof(int) * (l1 + 1));
        (*M_score)[i] = malloc(sizeof(int) * (l1 + 1));
    }

    for (int row = 1; row < l2 + 1; row++){
        X2[row][0] = -params.GO + (-params.GE * (row - 1));
        Y2[row][0] = -1000;
        MM[row][0] = -1000;
        B0[row][0] = 1;
        B1[row][0] = 1;
        B2[row][0] = 1;
    }
    for (int col = 1; col < l1 + 1; col++){
        X2[0][col] = -1000;
        Y2[0][col] = -params.GO + (-params.GE * (col - 1));
        MM[0][col] = -1000;
        B0[0][col] = 2;
        B1[0][col] = 2;
        B2[0][col] = 2;
        
    }

    MM[0][0] = 0;
    X2[0][0] = 0;
    Y2[0][0] = 0;

    int s;
    
    for(int row = 1; row < l2 + 1; row++){
        for(int col = 1; col < l1 + 1; col++){
            
            s = mmMatrix(params.mmflag, s1[col-1],s2[row-1]);
          
            MM[row][col] = max_two_number(max_two_number(MM[row-1][col-1] + s, X2[row-1][col-1] + s), Y2[row-1][col-1] + s);
            X2[row][col] = max_two_number(max_two_number(MM[row-1][col] - params.GO, X2[row-1][col] - params.GE), Y2[row-1][col]-params.GO);
            Y2[row][col] = max_two_number(max_two_number(MM[row][col-1] - params.GO, Y2[row][col-1] - params.GE), X2[row][col-1]-params.GO);
            (*M_score)[row][col] = MM[row][col];
            //printf("%d row:%d col:%d * ", MM[row][col], row, col);

            if (MM[row][col] == MM[row-1][col-1]+s)
                B0[row][col] = 0;
            else if (MM[row][col] == X2[row-1][col-1]+s)
                B0[row][col] = 1;
            else 
                B0[row][col] = 2;

            if(MM[row][col] == MM[row-1][col-1]+s && MM[row][col] == X2[row-1][col-1]+s)
                B0[row][col] = 1;

            if(MM[row][col] == MM[row-1][col-1]+s && MM[row][col] == Y2[row-1][col-1]+s)
                B0[row][col] = 2;

            if(MM[row][col] == X2[row-1][col-1]+s && MM[row][col] == Y2[row-1][col-1]+s)
                B0[row][col] = 3;
            
            if(X2[row][col] == MM[row-1][col]-params.GO)
                B1[row][col] = 0;
            else if(X2[row][col] == X2[row-1][col]-params.GE)
                B1[row][col] = 1;
            else if(X2[row][col] == Y2[row-1][col]-params.GO)
                B1[row][col] = 2;

            if(Y2[row][col] == MM[row][col-1]-params.GO)
                B2[row][col] = 0;
            else if(Y2[row][col] == X2[row][col-1]-params.GO)
                B2[row][col] = 1;
            else if(Y2[row][col] == Y2[row][col-1]-params.GE)
                B2[row][col] = 2;
        }
    }

    
    int R = l2;
    int C = l1;

    int max = max_two_number(max_two_number(MM[R][C], X2[R][C]), Y2[R][C]);

    if (max == MM[R][C])
        B0_func(B0, B1, B2, Aligned_pair, s1, s2, C, R, ct1, ct2, final_length);
    else if(max == X2[R][C])
        B1_func(B0, B1, B2, Aligned_pair, s1, s2, C, R, ct1, ct2, final_length);
    else if(max == Y2[R][C])
        B2_func(B0, B1, B2, Aligned_pair, s1, s2, C, R, ct1, ct2, final_length);


    //*M_score = MM;
    *X_score = X2;
    *Y_score = Y2;
    //printf("%d\n\n", (*M_score)[3][3]);
    

}

int Alignment_score_2D(int* seq_1, int* seq_2, int l, AlignmentParams params){

    int score = 0;
    int gflag1 = 1;
    int gflag2 = 1;

    for(int i=0; i<l; i++){

        if(seq_1[i] != 0 && seq_2[i] != 0){
            score = score + mmMatrix( params.mmflag, seq_1[i], seq_2[i]);
            gflag1 = 1;
            gflag2 = 1;
        }else if(seq_1[i] == 0 && seq_2[i] != 0){
            gflag2 = 1;
            if(gflag1 == 1){
                gflag1 = 0;
                score = score - params.GO;
            }
            else
                score = score - params.GE;
        } else if(seq_1[i] != 0 && seq_2[i] == 0){
            gflag1 = 1;
            if(gflag2 == 1){
                gflag2 = 0;
                score = score - params.GO;
            }
            else
                score = score - params.GE;
        }
        else{
            gflag1 = 0;
            gflag2 = 0;
        }
    }

    return score;
}

void Alignment_scoring_matrix(FILE* output, int* seq1, int* seq2, int* rev_seq1, int* rev_seq2, 
    int seq1_length, int seq2_length, int*** Aligned_pair, int*** score_matrix, AlignmentParams params, int* score, int* alignment_2D_l){

    int** f_MM = NULL;
    int** f_X2 = NULL;
    int** f_Y2 = NULL;
    int** f_B0 = NULL;
    int** r_MM = NULL;
    int** r_X2 = NULL;
    int** r_Y2 = NULL;
    
    int alignment_2D_length = 0;

    
    
    Alignment_2D(Aligned_pair, &f_MM, &f_X2, &f_Y2, seq1, seq2, seq1_length, seq2_length, &alignment_2D_length, params);

    //fprintf(output, "the length of 2D alignment : %d\n\n", alignment_2D_length);
    //fprintf(output, "sequence1\n");
    //write_sequence(output, alignment_2D_length, (*Aligned_pair)[0]);
    //fprintf(output, "sequence2\n");
    //write_sequence(output, alignment_2D_length, (*Aligned_pair)[1]);

    alignment_2D_length = 0;
    Alignment_2D(Aligned_pair, &r_MM, &r_X2, &r_Y2, rev_seq1, rev_seq2, seq1_length, seq2_length, &alignment_2D_length, params);

    //fprintf(output, "the length of 2D alignment : %d\n\n", alignment_2D_length);
    //fprintf(output, "sequence1\n");
    //write_sequence(output, alignment_2D_length, (*Aligned_pair)[0]);
    //fprintf(output, "sequence2\n");
    //write_sequence(output, alignment_2D_length, (*Aligned_pair)[1]);

    int c1 = 0;
    int c2 = 0;

    int a,b,c;
    int ip=0, jp=0;
    for(int i = 1; i<seq2_length ; i++){
        for(int j = 1; j<seq1_length ; j++){

            ip = seq2_length - (i);
            jp = seq1_length - (j);

            //a = max_two_number(max_two_number(f_MM[i][j]+r_MM[ip][jp], f_Y2[i][j] + r_Y2[ip][jp]+GO-GE), (f_X2[i][j] + r_X2[ip][jp]+GO-GE));
            //b = max_two_number(max_two_number(f_MM[i][j]+r_X2[ip][jp], f_MM[i][j] + r_Y2[ip][jp]), f_X2[i][j] + r_Y2[ip][jp]);
            //c = max_two_number(max_two_number(f_X2[i][j]+r_MM[ip][jp], f_Y2[i][j] +r_MM[ip][jp]), f_Y2[i][j] + r_X2[ip][jp]);
            //(*score_matrix)[i][j] = max_two_number(max_two_number(a,b),c);
            (*score_matrix)[i][j] = f_MM[i][j] + r_MM[ip][jp];
            //fprintf(output, "%d %d\n%d %d\n%d\n",i,j,f_MM[i][j], r_MM[ip][jp], (*score_matrix)[i][j]);

        }
    }
    for(int i = 1; i < seq1_length + 1; i ++){
        (*score_matrix)[0][i] = (*score_matrix)[1][1];
    }
    for(int i = 1; i < seq2_length + 1; i ++){
        (*score_matrix)[i][0] = (*score_matrix)[1][1];
    }
    *alignment_2D_l = alignment_2D_length;
    *score = Alignment_score_2D((*Aligned_pair)[0], (*Aligned_pair)[1], alignment_2D_length, params);
    //printf("%d\n***\n", *score);
    
}