# 3D Multiple Sequence Alignment with Affine Gap Penalties

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

This project implements a **three-way (3D) global sequence alignment** algorithm using dynamic programming and affine gap penalties. It is designed for protein sequences and can handle three sequences simultaneously, providing more accurate results than traditional pairwise methods.

## Features

- **3D Dynamic Programming:** Simultaneous alignment of three sequences.
- **Affine Gap Penalties:** Supports both gap opening and gap extension penalties.
- **Custom Scoring Matrices:** Includes BLOSUM62, BLOSUM30, and BLOSUM90.
- **Command-Line Parameters:** User-friendly prompts for penalties and matrix selection.
- **FASTA Input:** Reads protein sequences from a standard FASTA file.
- **Output:** Produces aligned sequences and alignment scores.

## Getting Started

### 1. Build

To build the program, you will need a C compiler (such as GCC).

**Compile using:**
```bash
gcc Readsequences.c dynamic_programming.c mmMatrix.c write_sequence.c Alignment_scoring_matrix.c  Triple.c -o Triple

### 2. Prepare Input
Prepare a file named input in the project directory, containing your protein sequences in FASTA format. Example:

>seq1
MKTAYIAKQRQISFVKSHFSRQDILDLWQ
>seq2
MKTAYIAKQ-KISFVKSHFSNQDILD--Q
>seq3
MKLAYIAKQRQISFVKSHFSNQDILD--Q

### 3. Run the Program
Run the executable:

./align3d
The program will prompt you for:

Gap opening penalty (GO)
Gap extension penalty (GE)
Scoring matrix choice (0: BLOSUM62, 1: BLOSUM30, 2: BLOSUM90)
### 4. Output
The aligned sequences and alignment scores will be written to alignment_output.txt.
A gapless version of the sequences will be written to a file named gapless.

File Structure

.
├── main.c
├── Readsequences.c
├── Readsequences.h
├── write_sequence.c
├── write_sequence.h
├── mmMatrix.c
├── mmMatrix.h
├── dynamic_programming.c
├── dynamic_programming.h
├── Alignment_params.c
├── Alignment_params.h
├── input
├── alignment_output.txt
├── gapless
├── README.md
└── LICENSE


### License

This project is licensed under the MIT License.

### Citation

If you use this code in your research, please cite:

Askari Rad, M. et al., "Three-Way Alignment Improves Multiple Sequence Alignment of Highly Diverged Sequences", Algorithms 2024.  https://www.mdpi.com/1999-4893/17/5/205
Contact

For questions, suggestions, or contributions, please open an issue or contact the maintainer:

Mahbubeh Askarirad
Email: mahbubeh.askarirad@gmail.com