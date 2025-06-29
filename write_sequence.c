#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "write_sequence.h"

static const int Amino_num[26] = {1,0,2,3,4,5,6,7,8,0,9,10,11,12,0,13,14,15,16,17,0,18,19,0,20,0};

void write_sequence(FILE* output, int l, int* seq){
    int count = 0;
    
    for(int i = 0; i < l; i++){
        if(seq[i] == 0) fprintf(output, "-");
        else {
            for(int j = 0; j < 25; j++){
                if(Amino_num[j] == seq[i]){
                    fprintf(output, "%c", j + 65);
                    break;
                }
            }
        }
        count++;
        if (count % 50 == 0){
            fprintf(output, "\n");
        }
    }
    if (count % 50 != 0) {
        fprintf(output, "\n");
    }
}

