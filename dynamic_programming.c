/**
 * Frees a 3D dynamically-allocated array of dimensions [dim1][dim2][...]
 * array: Pointer to the 3D array.
 * dim1: First dimension size.
 * dim2: Second dimension size.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mmMatrix.h"
#include "dynamic_programming.h"
#include "Alignment_params.h"

int*** M;
int*** XY;
int*** XZ;
int*** YZ;
int*** X;
int*** Y;
int*** Z;
int*** BM;
int*** BXY;
int*** BXZ;
int*** BYZ;
int*** BX;
int*** BY;
int*** BZ;

void free_3D(int*** array, int dim1, int dim2) {
    for (int i = 0; i < dim1; i++) {
        for (int j = 0; j < dim2; j++) {
            free(array[i][j]);
        }
        free(array[i]);
    }
    free(array);
}

void max_function(const int* nums, int length, int* max_func_out ){
    int max = nums[0];
    int idx = 0;
    for (int ctr = 1; ctr < length; ctr++ ){
        if (max < nums[ctr]){
            max = nums[ctr];
            idx = ctr;
        }    
    }
    max_func_out[0] = max;
    max_func_out[1] = idx;  
    
}

int BM_func(int* s1, int* s2, int* s3, int x, int y, int z, int*** AlignedTriplet, int* alignment_length);
int BXY_func(int* s1, int* s2, int* s3, int x, int y, int z, int*** AlignedTriplet, int* alignment_length);
int BXZ_func(int* s1, int* s2, int* s3, int x, int y, int z, int*** AlignedTriplet, int* alignment_length);
int BYZ_func(int* s1, int* s2, int* s3, int x, int y, int z, int*** AlignedTriplet, int* alignment_length);
int BX_func(int* s1, int* s2, int* s3, int x, int y, int z, int*** AlignedTriplet, int* alignment_length);
int BY_func(int* s1, int* s2, int* s3, int x, int y, int z, int*** AlignedTriplet, int* alignment_length);
int BZ_func(int* s1, int* s2, int* s3, int x, int y, int z, int*** AlignedTriplet, int* alignment_length);

int BM_func(int* s1, int* s2, int* s3, int x, int y, int z, int*** AlignedTriplet, int* alignment_length){
    if(x == 0 && y == 0 && z == 0 || x < 0 || y < 0 || z < 0)
        return 0;
    
    (*AlignedTriplet)[0][(*alignment_length)] = s1[x-1];
    (*AlignedTriplet)[1][(*alignment_length)] = s2[y-1];
    (*AlignedTriplet)[2][(*alignment_length)] = s3[z-1];
    
    //printf("M: %d %d %d\n%d %d\n",x , y, x, M[x][y][z], BM[x][y][z]);
    (*alignment_length)++;

    switch (BM[x][y][z])
    {
    case 0 :
        BM_func(s1, s2, s3, x-1, y-1, z-1, AlignedTriplet, alignment_length);
        break;
    case 1 :
        BXY_func(s1, s2, s3, x-1, y-1, z-1, AlignedTriplet, alignment_length);
        break;
    case 2 :
        BXZ_func(s1, s2, s3, x-1, y-1, z-1, AlignedTriplet, alignment_length);
        break;
    case 3 :
        BYZ_func(s1, s2, s3, x-1, y-1, z-1, AlignedTriplet, alignment_length);
        break;
    case 4 :
        BX_func(s1, s2, s3, x-1, y-1, z-1, AlignedTriplet, alignment_length);
        break;
    case 5 :
        BY_func(s1, s2, s3, x-1, y-1, z-1, AlignedTriplet, alignment_length);
        break;
    case 6 :
        BZ_func(s1, s2, s3, x-1, y-1, z-1, AlignedTriplet, alignment_length);
        break;
    
    default:
        break;
    }
    return 0;
}

int BXY_func(int* s1, int* s2, int* s3, int x, int y, int z, int*** AlignedTriplet, int* alignment_length){
    
    if(x == 0 && y == 0 && z == 0 || x < 0 || y < 0 || z < 0)
        return 0;

    (*AlignedTriplet)[0][(*alignment_length)] = s1[x-1];
    (*AlignedTriplet)[1][(*alignment_length)] = s2[y-1];
    (*AlignedTriplet)[2][(*alignment_length)] = 0;
    //printf("%d %d %d * %d %d %d\n", x,y,z, (*AlignedTriplet)[0][(*ctr1)], (*AlignedTriplet)[1][(*ctr2)], (*AlignedTriplet)[2][(*ctr3)]);

    //printf("XY: %d %d %d\n%d %d\n",x , y, x, XY[x][y][z], BXY[x][y][z]);
    (*alignment_length)++;

    switch (BXY[x][y][z])
    {
    case 0 :
       
        BM_func(s1, s2, s3, x-1, y-1, z, AlignedTriplet, alignment_length);
        break;
    case 1 :
        BXY_func(s1, s2, s3, x-1, y-1, z, AlignedTriplet, alignment_length);
        break;
    case 2 :
        BXZ_func(s1, s2, s3, x-1, y-1, z, AlignedTriplet, alignment_length);
        break;
    case 3 :
        BYZ_func(s1, s2, s3, x-1, y-1, z, AlignedTriplet, alignment_length);
        break;
    case 4 :
        BX_func(s1, s2, s3, x-1, y-1, z, AlignedTriplet, alignment_length);

        break;
    case 5 :
        BY_func(s1, s2, s3, x-1, y-1, z, AlignedTriplet, alignment_length);
        break;
    case 6 :
        BZ_func(s1, s2, s3, x-1, y-1, z, AlignedTriplet, alignment_length);
        break;
    
    default:
        break;
    }
    
    return 0;
}

int BXZ_func(int* s1, int* s2, int* s3, int x, int y, int z, int*** AlignedTriplet, int* alignment_length){
    
    if(x == 0 && y == 0 && z == 0 || x < 0 || y < 0 || z < 0)
        return 0;

    (*AlignedTriplet)[0][(*alignment_length)] = s1[x-1];
    (*AlignedTriplet)[1][(*alignment_length)] = 0;
    (*AlignedTriplet)[2][(*alignment_length)] = s3[z-1];
    //printf("%d %d %d * %d %d %d\n", x,y,z, (*AlignedTriplet)[0][(*ctr1)], (*AlignedTriplet)[1][(*ctr2)], (*AlignedTriplet)[2][(*ctr3)]);

    //printf("%d %d %d\n%d %d\n",x , y, x, XZ[x][y][z], BXZ[x][y][z]);
    (*alignment_length)++;

    switch (BXZ[x][y][z])
    {
    case 0 :
       
        BM_func(s1, s2, s3, x-1, y, z-1, AlignedTriplet, alignment_length);
        break;
    case 1 :
        BXY_func(s1, s2, s3, x-1, y, z-1, AlignedTriplet, alignment_length);
        break;
    case 2 :
        BXZ_func(s1, s2, s3, x-1, y, z-1, AlignedTriplet, alignment_length);
        break;
    case 3 :
        BYZ_func(s1, s2, s3, x-1, y, z-1, AlignedTriplet, alignment_length);
        break;
    case 4 :
        BX_func(s1, s2, s3, x-1, y, z-1, AlignedTriplet, alignment_length);
        break;
    case 5 :
        BY_func(s1, s2, s3, x-1, y, z-1, AlignedTriplet, alignment_length);
        break;
    case 6 :
        BZ_func(s1, s2, s3, x-1, y, z-1, AlignedTriplet, alignment_length);
        break;
    
    default:
        break;
    }
    
    return 0;
}

int BYZ_func(int* s1, int* s2, int* s3, int x, int y, int z, int*** AlignedTriplet, int* alignment_length){
    if(x == 0 && y == 0 && z == 0 || x < 0 || y < 0 || z < 0)
        return 0;

    (*AlignedTriplet)[0][(*alignment_length)] = 0;
    (*AlignedTriplet)[1][(*alignment_length)] = s2[y-1];
    (*AlignedTriplet)[2][(*alignment_length)] = s3[z-1];
    //printf("%d %d %d * %d %d %d\n", x,y,z, (*AlignedTriplet)[0][(*ctr1)], (*AlignedTriplet)[1][(*ctr2)], (*AlignedTriplet)[2][(*ctr3)]);

    //printf("YZ: %d %d %d\n%d %d\n",x , y, x, YZ[x][y][z], BYZ[x][y][z]);
    (*alignment_length)++;

    switch (BYZ[x][y][z])
    {
    case 0 :
       
        BM_func(s1, s2, s3, x, y-1, z-1, AlignedTriplet, alignment_length);
        break;
    case 1 :
        BXY_func(s1, s2, s3, x, y-1, z-1, AlignedTriplet, alignment_length);
        break;
    case 2 :
        BXZ_func(s1, s2, s3, x, y-1, z-1, AlignedTriplet, alignment_length);
        break;
    case 3 :
        BYZ_func(s1, s2, s3, x, y-1, z-1, AlignedTriplet, alignment_length);
        break;
    case 4 :
        BX_func(s1, s2, s3, x, y-1, z-1, AlignedTriplet, alignment_length);
        break;
    case 5 :
        BY_func(s1, s2, s3, x, y-1, z-1, AlignedTriplet, alignment_length);
        break;
    case 6 :
        BZ_func(s1, s2, s3, x, y-1, z-1, AlignedTriplet, alignment_length);
        break;
    
    default:
        break;
    }
    return 0;
}

int BX_func(int* s1, int* s2, int* s3, int x, int y, int z, int*** AlignedTriplet, int* alignment_length){
    if(x == 0 && y == 0 && z == 0 || x < 0 || y < 0 || z < 0)
        return 0;

    (*AlignedTriplet)[0][(*alignment_length)] = s1[x-1];
    (*AlignedTriplet)[1][(*alignment_length)] = 0;
    (*AlignedTriplet)[2][(*alignment_length)] = 0;
    //printf("%d %d %d * %d %d %d\n", x,y,z, (*AlignedTriplet)[0][(*ctr1)], (*AlignedTriplet)[1][(*ctr2)], (*AlignedTriplet)[2][(*ctr3)]);
    
    (*alignment_length)++;
    //printf("X: %d %d %d\n%d %d\n",x , y, x, X[x][y][z], BX[x][y][z]);

    switch (BX[x][y][z])
    {
    case 0 :
        BM_func(s1, s2, s3, x-1, y, z, AlignedTriplet, alignment_length);
        break;
    case 1 :
        BXY_func(s1, s2, s3, x-1, y, z, AlignedTriplet, alignment_length);
        break;
    case 2 :
        BXZ_func(s1, s2, s3, x-1, y, z, AlignedTriplet, alignment_length);
        break;
    case 3 :
        BYZ_func(s1, s2, s3, x-1, y, z, AlignedTriplet, alignment_length);
        break;
    case 4 :
        BX_func(s1, s2, s3, x-1, y, z, AlignedTriplet, alignment_length);
        break;
    case 5 :
        BY_func(s1, s2, s3, x-1, y, z, AlignedTriplet, alignment_length);
        break;
    case 6 :
        BZ_func(s1, s2, s3, x-1, y, z, AlignedTriplet, alignment_length);
        break;
    
    default:
        break;
    }
    return 0;
}

int BY_func(int* s1, int* s2, int* s3, int x, int y, int z, int*** AlignedTriplet, int* alignment_length){
    if(x == 0 && y == 0 && z == 0 || x < 0 || y < 0 || z < 0)
        return 0;

    (*AlignedTriplet)[0][(*alignment_length)] = 0;
    (*AlignedTriplet)[1][(*alignment_length)] = s2[y-1];
    (*AlignedTriplet)[2][(*alignment_length)] = 0;
    //printf("%d %d %d * %d %d %d\n", x,y,z, (*AlignedTriplet)[0][(*ctr1)], (*AlignedTriplet)[1][(*ctr2)], (*AlignedTriplet)[2][(*ctr3)]);
    //printf("Y: %d %d %d\n%d %d\n",x , y, x, Y[x][y][z], BY[x][y][z]);
    (*alignment_length)++;

    switch (BY[x][y][z])
    {
    case 0 :
       
        BM_func(s1, s2, s3, x, y-1, z, AlignedTriplet, alignment_length);
        break;
    case 1 :
        BXY_func(s1, s2, s3, x, y-1, z, AlignedTriplet, alignment_length);
        break;
    case 2 :
        BXZ_func(s1, s2, s3, x, y-1, z, AlignedTriplet, alignment_length);
        break;
    case 3 :
        BYZ_func(s1, s2, s3, x, y-1, z, AlignedTriplet, alignment_length);
        break;
    case 4 :
        BX_func(s1, s2, s3, x, y-1, z, AlignedTriplet, alignment_length);
        break;
    case 5 :
        BY_func(s1, s2, s3, x, y-1, z, AlignedTriplet, alignment_length);
        break;
    case 6 :
        BZ_func(s1, s2, s3, x, y-1, z, AlignedTriplet, alignment_length);
        break;
    
    default:
        break;
    }
    return 0;
}

int BZ_func(int* s1, int* s2, int* s3, int x, int y, int z, int*** AlignedTriplet, int* alignment_length){
if(x == 0 && y == 0 && z == 0 || x < 0 || y < 0 || z < 0)
        return 0;

    (*AlignedTriplet)[0][(*alignment_length)] = 0;
    (*AlignedTriplet)[1][(*alignment_length)] = 0;
    (*AlignedTriplet)[2][(*alignment_length)] = s3[z-1];
    //printf("%d %d %d * %d %d %d\n", x,y,z, (*AlignedTriplet)[0][(*ctr1)], (*AlignedTriplet)[1][(*ctr2)], (*AlignedTriplet)[2][(*ctr3)]);

    //printf("Z: %d %d %d\n%d %d\n",x , y, x, Z[x][y][z], BZ[x][y][z]);
    (*alignment_length)++;

    switch (BZ[x][y][z])
    {
    case 0 :
        BM_func(s1, s2, s3, x, y, z-1, AlignedTriplet, alignment_length);
        break;
    case 1 :
        BXY_func(s1, s2, s3, x, y, z-1, AlignedTriplet, alignment_length);
        break;
    case 2 :
        BXZ_func(s1, s2, s3, x, y, z-1, AlignedTriplet, alignment_length);
        break;
    case 3 :
        BYZ_func(s1, s2, s3, x, y, z-1, AlignedTriplet, alignment_length);
        break;
    case 4 :
        BX_func(s1, s2, s3, x, y, z-1, AlignedTriplet, alignment_length);
        break;
    case 5 :
        BY_func(s1, s2, s3, x, y, z-1, AlignedTriplet, alignment_length);
        break;
    case 6 :
        BZ_func(s1, s2, s3, x, y, z-1, AlignedTriplet, alignment_length);
        break;
    
    default:
        break;
    }
    return 0;
}

int*** allocate_3D(int dim1, int dim2, int dim3) {
    int*** array = malloc(dim1 * sizeof(int**));
    if (!array) {
        fprintf(stderr, "Error: malloc failed for 3D array (dim1)\n");
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < dim1; i++) {
        array[i] = malloc(dim2 * sizeof(int*));
        if (!array[i]) {
            fprintf(stderr, "Error: malloc failed for 3D array (dim2)\n");
            // Free already allocated memory before exit
            for (int x = 0; x < i; x++) free(array[x]);
            free(array);
            exit(EXIT_FAILURE);
        }
        for (int j = 0; j < dim2; j++) {
            array[i][j] = malloc(dim3 * sizeof(int));
            if (!array[i][j]) {
                fprintf(stderr, "Error: malloc failed for 3D array (dim3)\n");
                // Free already allocated memory before exit
                for (int y = 0; y < j; y++) free(array[i][y]);
                for (int x = 0; x < i; x++) {
                    for (int y = 0; y < dim2; y++) free(array[x][y]);
                    free(array[x]);
                }
                free(array[i]);
                free(array);
                exit(EXIT_FAILURE);
            }
        }
    }
    return array;
}

void dynamic_programming(int* s1, int* s2, int* s3, int seq1_length, int seq2_length, int seq3_length, int*** AlignedTriplet, AlignmentParams params, int* triplet_length){
    
    int max_func_out[2];
    int alignment_length = 0;
    int values[7];
    int score = 0;
    int score12, score13, score23;
    *AlignedTriplet = malloc(sizeof(int*) * 3);
    for (int k=0; k<3; k++){
        (*AlignedTriplet)[k] = malloc(sizeof(int) * (seq1_length + seq2_length + seq3_length));
    }

    M = allocate_3D(seq1_length + 1, seq2_length + 1, seq3_length + 1);
    XY = allocate_3D(seq1_length + 1, seq2_length + 1, seq3_length + 1);
    XZ = allocate_3D(seq1_length + 1, seq2_length + 1, seq3_length + 1);
    YZ = allocate_3D(seq1_length + 1, seq2_length + 1, seq3_length + 1);
    X = allocate_3D(seq1_length + 1, seq2_length + 1, seq3_length + 1);
    Y = allocate_3D(seq1_length + 1, seq2_length + 1, seq3_length + 1);
    Z = allocate_3D(seq1_length + 1, seq2_length + 1, seq3_length + 1);
    BM = allocate_3D(seq1_length + 1, seq2_length + 1, seq3_length + 1);
    BXY = allocate_3D(seq1_length + 1, seq2_length + 1, seq3_length + 1);
    BXZ = allocate_3D(seq1_length + 1, seq2_length + 1, seq3_length + 1);
    BYZ = allocate_3D(seq1_length + 1, seq2_length + 1, seq3_length + 1);
    BX = allocate_3D(seq1_length + 1, seq2_length + 1, seq3_length + 1);
    BY = allocate_3D(seq1_length + 1, seq2_length + 1, seq3_length + 1);
    BZ = allocate_3D(seq1_length + 1, seq2_length + 1, seq3_length + 1);
    
    /*
    for(int i=0; i<seq1_length + 1; i++){printf("ERROR3\n");
        for(int j =0; j< seq2_length + 1; j++){
            for(int k=0; k<seq3_length + 1; k++){ 

                M[i][j][k] = -100000;
                XY[i][j][k] = -100000;
                XZ[i][j][k] = -100000;
                YZ[i][j][k] = -100000;
                X[i][j][k] = -100000;
                Y[i][j][k] = -100000;
                Z[i][j][k] = -100000;
                BXY[i][j][k] = -100000;
                BXZ[i][j][k] = -100000;
                BYZ[i][j][k] = -100000;
                BX[i][j][k] = -100000;
                BY[i][j][k] = -100000;
                BZ[i][j][k] = -100000;
                BM[i][j][k] = -100000;
            }
        }
    }
*/
    
    
    int cc = 0;
    int cube_capacity = (seq1_length+1) * (seq2_length+1) * (seq3_length+1);
    printf("the total number of searching cells: %d\n", cube_capacity);
    printf("Dynamic programing started ...\n");
for(int i=0; i<seq1_length + 1; i++){
        for(int j =0; j< seq2_length + 1; j++){
            for(int k=0; k<seq3_length + 1; k++){

                M[i][j][k] = -100000;
                XY[i][j][k] = -100000;
                XZ[i][j][k] = -100000;
                YZ[i][j][k] = -100000;
                X[i][j][k] = -100000;
                Y[i][j][k] = -100000;
                Z[i][j][k] = -100000;
                BXY[i][j][k] = -100000;
                BXZ[i][j][k] = -100000;
                BYZ[i][j][k] = -100000;
                BX[i][j][k] = -100000;
                BY[i][j][k] = -100000;
                BZ[i][j][k] = -100000;
                BM[i][j][k] = -100000;
            }
        }
    }

   
     M[0][0][0] = 0;
    XY[0][0][0] = 0;
    XZ[0][0][0] = 0;
    YZ[0][0][0] = 0;
    X[0][0][0] = 0;
    Y[0][0][0] = 0;
    Z[0][0][0] = 0;

    for(int z = 0; z < seq3_length + 1; z++){
        for(int y = 0; y < seq2_length + 1; y++){
            for(int x = 0; x < seq1_length + 1; x++){
                cc++;
                if(x == 0 || y == 0 || z == 0){
                // NAN
                }
                if(x == 1 && y == 1 && z == 1){ // first value // there sequences aligned with their first base
                    score12 = mmMatrix(params.mmflag, s1[x-1],s2[y-1]);
                    score13 = mmMatrix(params.mmflag, s1[x-1],s3[z-1]);
                    score23 = mmMatrix(params.mmflag, s2[y-1],s3[z-1]);
                    M[1][1][1] = 0 + score12 + score13 + score23;
                    BM[1][1][1] = 0;
                }
                
                if(x > 1 && y > 1 && z > 1 || x == 1 && y > 1 && z > 1 || x > 1 && y == 1 && z > 1 ||
                x > 1 && y > 1 && z == 1 || x == 1 && y == 1 && z > 1 || x == 1 && y > 1 && z == 1 ||
                x > 1 && y == 1 && z == 1){
                    score12 = mmMatrix(params.mmflag, s1[x-1],s2[y-1]);
                    score13 = mmMatrix(params.mmflag, s1[x-1],s3[z-1]);
                    score23 = mmMatrix(params.mmflag, s2[y-1],s3[z-1]);
                    score = score12 + score13+ score23; 
                    values[0] = M[x-1][y-1][z-1] + score;
                    values[1] = XY[x-1][y-1][z-1] + score;
                    values[2] = XZ[x-1][y-1][z-1] + score;
                    values[3] = YZ[x-1][y-1][z-1] + score;
                    values[4] = X[x-1][y-1][z-1] + score;
                    values[5] = Y[x-1][y-1][z-1] + score;
                    values[6] = Z[x-1][y-1][z-1] + score;

                    max_function(values, 7, max_func_out);
                    M[x][y][z] = max_func_out[0];
                    BM[x][y][z] = max_func_out[1]; 
                }
                // XY[x][y][z] *********
                if(x == 0 || y == 0){
                    // NAN
                }
                if(x==1 && y==1 && z==0){ // the first value when S1 and s2 liparams.GEnd with thier first base and s3 has gap
                    score12 = mmMatrix(params.mmflag, s1[x-1],s2[y-1]);
                    XY[1][1][0] = 0 - 2*params.GO + score12;
                    BXY[1][1][0] = 0;
                }
                if(x > 1 &&  y == 1 || x == 1 && y >1 || x > 1 && y > 1){
                    score12 = mmMatrix(params.mmflag, s1[x-1],s2[y-1]);
                    values[0] = M[x-1][y-1][z] - 2*params.GO + score12;
                    values[1] = XY[x-1][y-1][z] - 2*params.GE + score12;
                    values[2] = XZ[x-1][y-1][z] - 2*params.GO + score12;
                    values[3] = YZ[x-1][y-1][z] - 2*params.GO + score12;
                    values[4] = X[x-1][y-1][z] - 2*params.GE + score12;
                    values[5] = Y[x-1][y-1][z] - 2*params.GE + score12;
                    values[6] = Z[x-1][y-1][z] - 2*params.GO + score12;
                    
                    max_function(values, 7, max_func_out);
                    XY[x][y][z] = max_func_out[0];
                    BXY[x][y][z] = max_func_out[1];
                }
                // XZ
                if(x == 0 || z == 0){
                    // NAN
                }
                if(x==1 && y==0 && z == 1){// the first value 
                    score13 = mmMatrix(params.mmflag, s1[x-1],s3[z-1]);
                    XZ[1][0][1]= 0 - 2*params.GO + score13;
                    BXZ[1][0][1] = 0;
                }
                if(x > 1 && z == 1 || x == 1 && z > 1 || x > 1 && z > 1){
                    score13 = mmMatrix(params.mmflag, s1[x-1],s3[z-1]);
                    values[0] = M[x-1][y][z-1] - 2*params.GO + score13;
                    values[1] = XY[x-1][y][z-1] - 2*params.GO + score13;
                    values[2] = XZ[x-1][y][z-1] - 2*params.GE + score13;
                    values[3] = YZ[x-1][y][z-1] - 2*params.GO + score13;
                    values[4] = X[x-1][y][z-1] - 2*params.GE + score13;
                    values[5] = Y[x-1][y][z-1] - 2*params.GO + score13;
                    values[6] = Z[x-1][y][z-1] - 2*params.GE + score13;
                    max_function(values, 7, max_func_out);
                    XZ[x][y][z] = max_func_out[0];
                    BXZ[x][y][z] = max_func_out[1];

                }
                // YZ
                if(y == 0 || z == 0){
                    // NAN
                }
                if(x==0 && y==1 && z==1){ // first value
                    score23 = mmMatrix(params.mmflag, s2[y-1],s3[z-1]);
                    YZ[0][1][1] = 0 - 2*params.GO + score23;
                    BYZ[0][1][1] = 0;
                }
                if(y > 1 && z == 1|| y == 1 && z > 1 || y > 1 && z > 1){
                    score23 = mmMatrix(params.mmflag, s2[y-1],s3[z-1]);
                    values[0] = M[x][y-1][z-1] - 2*params.GO + score23;
                    values[1] = XY[x][y-1][z-1] - 2*params.GO + score23;
                    values[2] = XZ[x][y-1][z-1] - 2*params.GO + score23;
                    values[3] = YZ[x][y-1][z-1] - 2*params.GE + score23;
                    values[4] = X[x][y-1][z-1] - 2*params.GO + score23;
                    values[5] = Y[x][y-1][z-1] - 2*params.GE + score23;
                    values[6] = Z[x][y-1][z-1] - 2*params.GE + score23;
                    max_function(values, 7, max_func_out);
                    YZ[x][y][z] = max_func_out[0];
                    BYZ[x][y][z] = max_func_out[1];
                }
                // X
                if(x == 0){
                    //NAN
                }
                if(x==1 && y==0 && z == 0){ //first value
                    X[1][0][0] = 0 - 2*params.GO;
                    BX[1][0][0] = 0;
                }
                if(x > 1){
                    values[0] = M[x-1][y][z] - 2 * params.GO;
                    values[1] = XY[x-1][y][z] - params.GO - params.GE;
                    values[2] = XZ[x-1][y][z] - params.GO - params.GE;
                    values[3] = YZ[x-1][y][z] - 2 * params.GO;
                    values[4] = X[x-1][y][z] - 2 * params.GE;
                    values[5] = Y[x-1][y][z] - params.GO - params.GE;
                    values[6] = Z[x-1][y][z] - params.GO - params.GE;
                    max_function(values, 7, max_func_out);
                    X[x][y][z] = max_func_out[0];
                    BX[x][y][z] = max_func_out[1];
                }
                // Y  ************
                if(y == 0){
                    //Nan
                }
                if(x == 0  && y==1 && z == 0){ // first value 
                    
                    Y[0][1][0] = 0 - 2 * params.GO;
                    BY[0][1][0] = 0;
                }
                if(y > 1){
                    values[0] = M[x][y-1][z] - 2 * params.GO;
                    values[1] = XY[x][y-1][z] - params.GO - params.GE;
                    values[2] = XZ[x][y-1][z] - 2 * params.GO;
                    values[3] = YZ[x][y-1][z] - params.GO - params.GE;
                    values[4] = X[x][y-1][z] - params.GO - params.GE;
                    values[5] = Y[x][y-1][z] - 2 * params.GE;
                    values[6] = Z[x][y-1][z] - params.GO - params.GE;
                    max_function(values, 7, max_func_out);
                    Y[x][y][z] = max_func_out[0];
                    BY[x][y][z] = max_func_out[1];
                }
                // Z
                if(z == 0){
                    //NAN
                }
                if(x == 0 && y == 0 && z==1) { // first value
                    Z[0][0][1] =  0 - 2 * params.GO;
                    BZ[0][0][1] = 0;
                }   
                if(z > 1){
                    values[0] = M[x][y][z-1] - 2 * params.GO;
                    values[1] = XY[x][y][z-1] - 2 * params.GO;
                    values[2] = XZ[x][y][z-1] - params.GO - params.GE;
                    values[3] = YZ[x][y][z-1] - params.GO - params.GE;
                    values[4] = X[x][y][z-1] - params.GO - params.GE;
                    values[5] = Y[x][y][z-1] - params.GO - params.GE;
                    values[6] = Z[x][y][z-1] - 2 * params.GE;
                    max_function(values, 7, max_func_out);
                    Z[x][y][z] = max_func_out[0];
                    BZ[x][y][z] = max_func_out[1];
                }
            }
        }
    }
    
    printf("scoring matrices and traceback matrices been calculated...!\n");

    int x = seq1_length;
    int y = seq2_length;
    int z = seq3_length;

    values[0] = M[x][y][z];
    values[1] = XY[x][y][z];
    values[2] = XZ[x][y][z];
    values[3] = YZ[x][y][z];
    values[4] = X[x][y][z];
    values[5] = Y[x][y][z];
    values[6] = Z[x][y][z];
    
    
    max_function(values, 7, max_func_out);
    printf("%d %d\n", max_func_out[0], max_func_out[1]);
   
    printf("Traceback started ....\n");
    switch (max_func_out[1])
    { 
    case 0 :
        BM_func(s1, s2, s3, x, y, z, AlignedTriplet, &alignment_length);
        break;
    case 1 :
        BXY_func(s1, s2, s3, x, y, z, AlignedTriplet, &alignment_length);
        break;
    case 2 :
        BXZ_func(s1, s2, s3, x, y, z, AlignedTriplet, &alignment_length);
        break;
    case 3 :
        BYZ_func(s1, s2, s3, x, y, z, AlignedTriplet, &alignment_length);
        break;
    case 4 :
        BX_func(s1, s2, s3, x, y, z, AlignedTriplet, &alignment_length);
        break;
    case 5 :
        BY_func(s1, s2, s3, x, y, z, AlignedTriplet,  &alignment_length);
        break;
    case 6 :
        BZ_func(s1, s2, s3, x, y, z, AlignedTriplet, &alignment_length);
        break;
    
    default:
        break;
    }
    printf("Traceback is done!\n");

    *triplet_length = alignment_length;

    int dim1 = seq1_length + 1;
    int dim2 = seq2_length + 1;
    // Use these for all matrices:
    free_3D(M, dim1, dim2);
    free_3D(XY, dim1, dim2);
    free_3D(XZ, dim1, dim2);
    free_3D(YZ, dim1, dim2);
    free_3D(X, dim1, dim2);
    free_3D(Y, dim1, dim2);
    free_3D(Z, dim1, dim2);
    free_3D(BM, dim1, dim2);
    free_3D(BXY, dim1, dim2);
    free_3D(BXZ, dim1, dim2);
    free_3D(BYZ, dim1, dim2);
    free_3D(BX, dim1, dim2);
    free_3D(BY, dim1, dim2);
    free_3D(BZ, dim1, dim2);
    
}

