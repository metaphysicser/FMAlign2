#include <stddef.h>

typedef short int bool;


long long right_band(int n, int m, int i, int band);
long long left_band(int n, int m, int i, int band);
bool insideband(int x, int y, int band, int n, int m); 
double protein_score(char a, char b, double **mtx);
size_t query_place(int x, int y, int band, int n, int m);
int matrix_set(int x, int y, int band, double val, double *mtx, int n, int m);
double matrix_query(int x, int y, int band, double *mtx, int n, int m);
int matrix_update(int x, int y, int band, double val, double *mtx, int n, int m);
int matrix_set_INT(int x, int y, int band, int val, int *mtx, int n, int m);
int matrix_query_INT(int x, int y, int band, int *mtx, int n, int m);
