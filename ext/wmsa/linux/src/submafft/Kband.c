#include "Kband.h"
#include <stdio.h>
#include <stddef.h>

#define DEBUG 0

extern int nalphabets, amino[];
#define INF (1e+100)

long long __min__(long long x, long long y) { return x > y ? y : x; }
long long __max__(long long x, long long y) { return x > y ? x : y; }
long long __abs__(long long x)              { return x >= 0ll ? x : -x; }


long long right_band(int n, int m, int i, int band)
{
	return __min__(i + m - n + band, m - 1);
}

long long left_band(int n, int m, int i, int band)
{
	return __max__(i - band, 0);
}

bool insideband(int x, int y, int band, int n, int m)
{
    /* Judge the coordinate in the band */
	if(x < 0 || y < 0 || x >= n || y >= m) return 0;
	// abs(y - (x + (m - n) / 2)) <= band + (m - n) / 2
	// reporterr("\033[0;31mx = %d y = %d band = %d n = %d m = %d ? %d\033[0m\n", x, y, band, n, m, __abs__(y * 2 - (x * 2 + m - n)) <= band * 2 + m - n);
	return __abs__(y * 2 - (x * 2 + m - n)) <= band * 2 + m - n;
}

double protein_score(char a, char b, double **mtx)
{
	/* calcuate the score of c1 and c2, using BLOSUM62 */
	return mtx[a][b];
}

size_t query_place(int x, int y, int band, int n, int m)
{
	/* Query the place in the array 
	 * Now the complexity of the segment is O(1)
	 * see https://gitee.com/wym6912/grad_learn/blob/master/kband_mapping/mapping_rectangle_O1.py
	 */
	if(! insideband(x, y, band, n, m)) return -1;
	else 
	{
		if(x == 0) return (size_t)__abs__(y - x);
		size_t longest = right_band(n, m, n / 2, band) - left_band(n, m, n / 2, band) + 1, \
	           firstline = right_band(n, m, 0, band) - left_band(n, m, 0, band) + 1;
		if(firstline == longest)
		{
#if DEBUG
			fprintf(stderr, "A %lld\n", (longest * x + y - left_band(n, m, 0, band)));
#endif
			return (size_t)(longest * x + y - left_band(n, m, 0, band));
		}
		size_t ans = 0, thisline = right_band(n, m, x, band) - left_band(n, m, x, band) + 1;
		size_t beforeall, firstm, alllinesum, lastm, afterall;
		if(thisline == longest)
		{
			/* find the first line which has longest elements in the band */
			beforeall = longest - firstline + 1;
			ans += beforeall * (longest + firstline) / 2;
			firstm = longest - firstline; // get lines **before** this line which has longest elements
			alllinesum = x - 1 - firstm;
			ans += alllinesum * longest + y - left_band(n, m, x, band);
#if DEBUG
			fprintf(stderr, "B %llu\n", ans);
#endif			
			return ans;
		}
		if(thisline - firstline == x) // on [0, band]
		{
			ans = (thisline - 1 + firstline) * x / 2 + y - left_band(n, m, x, band);
#if DEBUG
			fprintf(stderr, "C %llu\n", ans);
#endif
			return ans;
		}
		beforeall = longest - firstline + 1;
		ans += beforeall * (longest + firstline) / 2;	
		firstm = longest - firstline;
		lastm = n - 1 - (longest - firstline);
		ans += longest * (lastm - firstm);
		afterall = longest - thisline - 1;
		ans += afterall * (longest + thisline) / 2 + y - left_band(n, m, x, band);
#if DEBUG
			fprintf(stderr, "D %llu\n", ans);
#endif
		return ans;
	}
	return -1;
}

int matrix_set(int x, int y, int band, double val, double *mtx, int n, int m)
{
	/* insert node by matrix */
	if(! insideband(x, y, band, n, m)) return -1;
	size_t q = query_place(x, y, band, n, m);
#if DEBUG
	fprintf(stderr, "(%d, %d) = %llu\n", x, y, q);
	fprintf(stderr, "%lf %lf\n", val, mtx[q]);
#endif
	if(q == -1) return -1;
	mtx[q] = val;
	return 0;
}

int matrix_update(int x, int y, int band, double val, double *mtx, int n, int m)
{
	/* insert node by matrix */
	if(! insideband(x, y, band, n, m)) return -1;
	size_t q = query_place(x, y, band, n, m);
#if DEBUG
	fprintf(stderr, "(%d, %d) = %llu; ", x, y, q);
	fprintf(stderr, "%lf %lf\n", val, mtx[q]);
#endif
	if(q == -1) return -1;
	mtx[q] += val;
#if DEBUG
	fprintf(stderr, "After update, mtx(%d, %d) = %lf\n", x, y, mtx[q]);
#endif
	return 0;
}

double matrix_query(int x, int y, int band, double *mtx, int n, int m)
{
	/* query node by matrix */
	size_t q = query_place(x, y, band, n, m);
	if(q == -1) return -INF;
	return mtx[q];
}

int matrix_set_INT(int x, int y, int band, int val, int *mtx, int n, int m)
{
	/* insert node by matrix
	 * type: int
	 */
	if(! insideband(x, y, band, n, m)) return -1;
	size_t q = query_place(x, y, band, n, m);
#if DEBUG
	fprintf(stderr, "(%d, %d) = %llu; ", x, y, q);
	fprintf(stderr, "%d %d\n", val, mtx[q]);
#endif
	if(q == -1) return -1;
	mtx[q] = val;
#if DEBUG
	fprintf(stderr, "After add, mtx(%d, %d) = %d\n", x, y, mtx[q]);
#endif
	return 0;
}

int matrix_query_INT(int x, int y, int band, int *mtx, int n, int m)
{
	/* query node by matrix 
	   type: int
	 */
	size_t q = query_place(x, y, band, n, m);
	if(q == -1) return -__INT_MAX__;
	return mtx[q];
}
