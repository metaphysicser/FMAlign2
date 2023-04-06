#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <sys/types.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <limits.h>
#include <stdarg.h>
#include "submafft/mtxutl.h"

#define N 5000000
#define B 256

int myatoi( char *in, char before ); // atoi + Expection on null string
double myatof( char *in, char before ); // atof + Expection on null string
void reporterr( const char *str, ... ); // Using stderr on error information
int myfgets(char s[], int l, FILE *fp); // fgets
int getnumlenandtype( FILE *fp, int seq_type ); // get the number and type of sequences, return the type of the sequences
int countKUorWA( FILE *fp ); // count the number of sequences
void readData( FILE *fp, char **name, int *nlen, char **seq, int njob ); // read fasta data with name and seq
void gapfilter_oneseq( char *aseq, char *seq ); // gap filter, write from seq to aseq
void writeData( FILE *fp, int locnjob, char **name, int *nlen, char **aseq ); // write data