#include <stdio.h> 
#include <unistd.h> 
#include <string.h>
#include <stdlib.h>

#define maxf  8192 // max file name length
#define maxn  500000 // max sequence length 

#define NOTKNOWNDOUBLE 1e+10
#define NOTKNOWNINT 1000007 

// filename
extern char *inputfile, *outputfile;

// sequence information
extern int seq_type; // 1 DNA, 2 Protein, 0 AUTO
extern int seq_num; // the number of sequences
extern int seq_maxlen; // the max length of sequences

// arguments on MAFFT
extern int fftWinSize, fftthreshold, BLOSUM, ppenalty, ppenalty_dist;
extern int alignband, threads, nmax_shift;
extern int printdebug;

// arguments on CD-HIT
extern double cdhitsim;
extern int maxmemory; // max memory for CD-HIT

// center sequence file name
extern char centerfile[];

// command stream
extern char cmdstr[], cmdstr2[], cmdstr3[];

// tmp directory
extern char tmpdir[], *readtmpdir;
extern int tmpinthisdir;

extern char orderprotein[], orderDNA[];
extern int BLOSUM62[20][20], trans[4][4];
