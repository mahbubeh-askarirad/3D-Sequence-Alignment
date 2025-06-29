/*
 * Readsequences.c
 * Reads protein sequences from a FASTA file into integer arrays using a custom encoding.
 * Also writes a gapless version of the input sequences to a file named "gapless".
 * Author: Mahbubeh Askarirad
 * Date: 2024-06-28
 * License: MIT
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Readsequences.h"

#define MAX_FILENAME_LENGTH 256
#define MAX_LINE_LEN 1000
#define MAX_SEQ_LENGTH 1000
#define MAX_SEQUENCES 3
#define AMINO_ALPHABET_SIZE 25
                                                // A   C,D,E,F,G,H,I, ,K,L ,M, N,   P, Q ,R ,S ,T ,  V, W
static const int Amino_num[AMINO_ALPHABET_SIZE] = {1,0,2,3,4,5,6,7,8,0,9,10,11,12,0,13,14,15,16,17,0,18,19,0,20};

// Decode amino code back to character
char decode_amino(int code) {
   for(int i = 0; i < 25; i++){
      if(Amino_num[i] == code) return 'A' + i;
   }
   return '?';
}
void Readsequences(FILE* output, const char* filename, int*** Seqs_out, int** lengths_out){ 
   FILE* fp = fopen(filename, "r");
   if (!fp){
      printf("failed to open input file");
      exit(EXIT_FAILURE);
   }
   printf("start reading the file ...\n");
   // Allocate temporary storage
   int* lengths = malloc(MAX_SEQUENCES * sizeof(int)); 
   if (!lengths) { fprintf(stderr, "Malloc failed for lengths\n"); exit(EXIT_FAILURE); }
   int** Seqs = malloc(MAX_SEQUENCES * sizeof(int*));
   if (Seqs == NULL) {
      fprintf(stderr, "Error: malloc failed for Seqs\n");
      free(lengths); // Free what you've already allocated
      exit(EXIT_FAILURE);
   }
    for (int i = 0; i < MAX_SEQUENCES; i++) {
        Seqs[i] = malloc(MAX_SEQ_LENGTH * sizeof(int));
        if (Seqs[i] == NULL) {
            fprintf(stderr, "Error: malloc failed for Seqs[%d]\n", i);
            // Free all previously allocated memory before exit
            for (int j = 0; j < i; j++) free(Seqs[j]);
            free(Seqs);
            free(lengths);
            exit(EXIT_FAILURE);
         }
    }
    

   char*   line = NULL;
   size_t  len = 0;
   ssize_t read;
   int row = -1, AA_counter = 0;

   while ((read = getline(&line, &len, fp)) != -1) {	
      if (line[read - 1] == '\n')  
         line[--read] = 0;	

      if (line[0] == '>'){
         if (row >= MAX_SEQUENCES - 1) {
            fprintf(stderr, "Error: more sequences than MAX_SEQUENCES (%d)\n", MAX_SEQUENCES);
            exit(EXIT_FAILURE);
         }
         AA_counter = 0;
         row++;
      }else{
         for(int counter = 0; line[counter] != '\0'; counter++){
            char c = line[counter];
            if(c >= 'A' && c <= 'Z'){
               int idx = c - 'A';
               if(Amino_num[idx] != 0){
                  if (AA_counter >= MAX_SEQ_LENGTH) {
                    fprintf(stderr, "Error: sequence too long (>%d)\n", MAX_SEQ_LENGTH);
                    exit(EXIT_FAILURE);
                  }
                  Seqs[row][AA_counter] = Amino_num[idx];
                  AA_counter++;
               }
               else{
                  fprintf(output,"Warning: unknown Amino Acid: %c !!", line[counter]);
               }
            }
            else if(c == '-'){
               // it is a gap;
            }
            else { 
               fprintf(output, "Warning: Invalid character '%c'\n", line[counter]);
                  continue;
               if(line[counter] != '-'){
               }
            }
         }
         lengths[row] = AA_counter;
      }
   }
   *Seqs_out = Seqs;
   *lengths_out = lengths;

   free(line);
   fclose(fp);

   // Write gapless output
   FILE* in = fopen(filename,"r");
   FILE* out = fopen("gapless","w");

   if (!in || !out) {
      printf("Error writing output");
      exit(EXIT_FAILURE);
   }

   char buf[MAX_LINE_LEN];
   int current = 0;
   while(fgets(buf, sizeof(buf), in)!=NULL){
      if (buf[0] == '>'){
         fputs(buf, out);
         for(int i = 0; i< lengths[current]; i++) {
            char amino = decode_amino(Seqs[current][i]);
            fputc(amino, out);
         }
         fputc('\n', out);
         current++;
      }
   }
   fclose(in);
   fclose(out);
}

// *