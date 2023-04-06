#include "mltaln.h"
#include "dp.h"
#include "mtxutl.h"
#include "Kband.h"

# define FFTNS
# define NOTPRINT
# define NOTPRINTLEN
# define DEBUG 0

# ifndef FFTNS //constants
const int nalphabets = 26;
const int nscoredalphabets = 20;
const int max_seq = 100000 + 7;
#endif

static const double inf = 1e100, eps = 1e-5;

double Kband(char *s1, char *s2, int len1, int len2, double penalty, int band, double **ad, char *r1, char *r2, int headgp, int tailgp)
{
	/*
		Every iteration of Kband.
		Return: Alignment score. -inf means this band can not finish this iteration.
	*/
	// fprintf(stderr, "In Kband, len1 = %d, len2 = %d, band = %d\n", len1, len2, band);
	//printf( "In Kband, len1 = %d, len2 = %d\n", len1, len2);
	int i, j, h, k;
	static TLS int len;
	static TLS char *ssss;
	/* 
	   insert node: matrix_set(key, val)
	   find node: map_t *data = matrix_query(&tree, key)
	   The complexity of this block is : time O(n), space O(n)
	*/
	static TLS int max_place_i = 0, moveval, tmpval_col, tmpval_row;
	static TLS double score = 0, tmpval, max_data_i = -inf;
	static TLS char *tmpr1, *tmpr2;
	static TLS int *move_vec, *max_place_j, *mpjp;
	static TLS double *val_vec, *max_data_j, *mdjp;
	move_vec = AllocateIntVecLarge((long long)((band << 1 | 1) + (len2 - len1)) * (len2 + 10) + 10);
	val_vec = AllocateDoubleVecLarge((long long)((band << 1 | 1) + (len2 - len1)) * (len2 + 10) + 10);
	max_place_j = AllocateIntVec(len2 + 10);
	max_data_j = AllocateDoubleVec(len2 + 10);

	/* Initial of the matrix */
	matrix_set(0, 0, band, protein_score(s1[0], s2[0], ad), val_vec, len1 + 1, len2 + 1);
	for(i = 1; i <= len1; ++ i) matrix_set(i, 0, band, protein_score(s1[i], s2[0], ad), val_vec, len1 + 1, len2 + 1);
	for(i = 1; i <= len2; ++ i) matrix_set(0, i, band, protein_score(s1[0], s2[i], ad), val_vec, len1 + 1, len2 + 1);
	
	if(headgp == 1)
	{
		for(i = 1; i <= len1; ++ i) matrix_update(i, 0, band, penalty, val_vec, len1 + 1, len2 + 1);
		for(i = 1; i <= len2; ++ i) matrix_update(0, i, band, penalty, val_vec, len1 + 1, len2 + 1);
	}

	/* Score Matrix and move_vec */
	max_data_i = -inf; max_place_i = 0; // Initial the max_data_i && max_place_i
	for(mdjp = max_data_j, mpjp = max_place_j, i = 0; i <= len2; ++ mdjp, ++ mpjp, ++ i) // Initial the max_data_j && max_place_j
	{
		*mdjp = matrix_query(0, i, band, val_vec, len1 + 1, len2 + 1);
		*mpjp = 0;
	}
	*max_data_j = 0;
	*max_place_j = 0;
	for(i = 1; i <= len1; ++ i)
	{
		max_data_i = matrix_query(i - 1, MAX(i - band, 0), band, val_vec, len1 + 1, len2 + 1); max_place_i = MAX(i - band, 0);
		if(insideband(i - 1, 0, band, len1 + 1, len2 + 1))
		{
			tmpval = matrix_query(i - 1, 0, band, val_vec, len1 + 1, len2 + 1);
			if(*max_data_j < tmpval)
			{
				*max_data_j = tmpval;
				*max_place_j = i - 1;
			}
		}	
		for(j = left_band(len1, len2, i, band); j <= right_band(len1 + 1, len2 + 1, i, band); ++ j)
		{
			if(1 <= j && j <= len2)
			{
				tmpval = matrix_query(i - 1, j - 1, band, val_vec, len1 + 1, len2 + 1);
				moveval = 0;
				/* Calcuate now the strategy of moving */
				if(insideband(i, j - 1, band, len1 + 1, len2 + 1)) 
				{
					if((tmpval_row = max_data_i + penalty) > tmpval)
					{
						tmpval = tmpval_row;
						moveval = -(j - max_place_i);
						// reporterr("Set by i, j - 1, j = %d, max_place_i = %d, val = %f\n", j, max_place_i, tmpval);
					}
					
				}
				if(insideband(i - 1, j, band, len1 + 1, len2 + 1))
				{
					if((tmpval_col = *(max_data_j + j) + penalty) > tmpval)
					{
						tmpval = tmpval_col;
						moveval = +(i - *(max_place_j + j));
						// reporterr("Set by i - 1, j, i = %d, *(max_place_j + j) = %d, val = %f\n", i, *(max_place_j + j), tmpval);
					}
				}
				if(i != len1 && j != len2) matrix_set(i, j, band, tmpval + protein_score(s1[i], s2[j], ad), val_vec, len1 + 1, len2 + 1);
				else matrix_set(i, j, band, tmpval, val_vec, len1 + 1, len2 + 1);
				matrix_set_INT(i, j, band, moveval, move_vec, len1 + 1, len2 + 1);

#if DEBUG
			reporterr("\n(%d, %d) = %d\nmax place i = %d, max data i = %f\nmax place j = ", i, j, matrix_query_INT(i, j, band, move_vec, len1 + 1, len2 + 1), max_place_i, max_data_i);
			for(k = 0; k <= len2; ++ k) reporterr("%6d ", *(max_place_j + k));
			reporterr("\nmax data j  = ");
			for(k = 0; k <= len2; ++ k) reporterr("%6.1f ", *(max_data_j + k));
			reporterr("\n\n");
#endif
				tmpval = matrix_query(i - 1, j - 1, band, val_vec, len1 + 1, len2 + 1);
				/* Calcuate the next */
				if(insideband(i, j - 1, band, len1 + 1, len2 + 1))
				{
					if(tmpval >= max_data_i)
					{
						max_data_i = tmpval;
						max_place_i = j - 1;
						// reporterr("Max_data_i: (%d, %d) -> %d %f\n", i, j - 1, max_place_i, max_data_i);
					}
				}
				if(insideband(i - 1, j, band, len1 + 1, len2 + 1))
				{
					if(tmpval >= *(max_data_j + j))
					{
						*(max_data_j + j) = tmpval;
						*(max_place_j + j) = i - 1;
						// reporterr("Max_data_j: (%d, %d) -> %d %f\n", i - 1, j, *(max_place_j + j), *(max_data_j + j));
					}
				}
			}
		}

		/* assure that the place appears in max_place_j is in the band */
		for(k = 0, mdjp = max_data_j, mpjp = max_place_j; k < i - band; ++ k, ++ mdjp, ++ mpjp) // clean the outside band date (only need to it?)
		{
			*mdjp = -inf + penalty;
			*mpjp = -__INT_MAX__;
		}
	}
	for(i = 0; i <= len1; ++ i) matrix_set_INT(i, 0, band, i + 1, move_vec, len1 + 1, len2 + 1);
	for(i = 0; i <= len2; ++ i) matrix_set_INT(0, i, band, -(i + 1), move_vec, len1 + 1, len2 + 1);

#if DEBUG
	reporterr("\n     ");
	for(j = 0; j < len2 + 1; ++ j, reporterr("%5d ", j - 1)); reporterr("\n%-5d", 0);
	for(i = 0; i < len1 + 1; ++ i, reporterr("\n%-5d", i))
		for(j = 0; j < len2 + 1; ++ j, reporterr(" "))
		{
			if(insideband(i, j, band, len1 + 1, len2 + 1))
			{
				reporterr("%5d", matrix_query_INT(i, j, band, move_vec, len1 + 1, len2 + 1));
			}
			else
			{
				reporterr("?????");
			}
		}
	reporterr("\n      ");
	for(j = 0; j < len2 + 1; ++ j, reporterr("%7d ", j - 1)); reporterr("\n%-6d", 0);
	for(i = 0; i < len1 + 1; ++ i, reporterr("\n%-6d", i))
		for(j = 0; j < len2 + 1; ++ j, reporterr(" "))
		{
			if(insideband(i, j, band, len1 + 1, len2 + 1))
			{
				reporterr("%7.1f", matrix_query(i, j, band, val_vec, len1 + 1, len2 + 1));
			}
			else
			{
				reporterr("???????");
			}
		}
	reporterr("\n");
#endif	

	/* Traceback */
	int limk = len1 + len2, nowi, nowj, val, tari, tarj;
	nowi = len1, nowj = len2;
	tmpr1 = r1 - 1, tmpr2 = r2 - 1;
	for(k = 0; k < limk; ++ k)
	{
		val = matrix_query_INT(nowi, nowj, band, move_vec, len1 + 1, len2 + 1);
		// calcuate the gaps
		if(val < 0) // column gap
		{
			tari = nowi - 1; tarj = nowj + val;
		}
		else if(val > 0) // row gap
		{
			tari = nowi - val; tarj = nowj - 1;
		}
		else // no gap
		{
			tari = nowi - 1; tarj = nowj - 1;
		}
#if DEBUG
		reporterr("(nowi, nowj) = (%d, %d) -> (tari, tarj) = (%d, %d), k = %d, limk = %d, val = %d\n", nowi, nowj, tari, tarj, k, limk, val);
#endif
		// insert it!
		h = nowi - tari;
		while(-- h > 0) // if the difference is lower or same than 1, than DO NOT insert it!
		{
			*++ tmpr1 = s1[tari + h];
			*++ tmpr2 = *newgapstr;  
			++ k;
		}
		h = nowj - tarj;
		while(-- h > 0)
		{
			*++ tmpr1 = *newgapstr;
			*++ tmpr2 = s2[tarj + h];
			++ k;
		}
		
		if(nowi <= 0 || nowj <= 0) break;
		*++ tmpr1 = s1[tari];
		*++ tmpr2 = s2[tarj];
		++ k;
		nowi = tari, nowj = tarj;
		// reporterr("r1 = %s, r2 = %s\n", r1, r2);
	}
	*++ tmpr1 = 0;
	*++ tmpr2 = 0;
	// reporterr("%s %s\n", r1, r2);
	score = matrix_query(len1, len2, band, val_vec, len1 + 1, len2 + 1);
	
	FreeIntVec(move_vec);
	FreeDoubleVec(val_vec);
	FreeIntVec(max_place_j);
	FreeDoubleVec(max_data_j);
	return score;
}

double G__align11( double **n_dynamicmtx, char **seq1, char **seq2, int alloclen, int headgp, int tailgp )
{
    //puts("In G__align11");
    /* 
        Kband Algorithm: align **two** sequences
        Return: Alignment Score. The sequence seq1 & seq2 must be aligned.
    */
    // Part -1: loop varibles
    int i, j;
# ifndef NOTPRINT
    // Part 0: Print some arguments
    printf("nalphabets = %d\n", nalphabets);
    for(i = 0; i < nalphabets; ++ i, puts(""))
        for(j = 0; j < nalphabets; ++ j, putchar(' '))
            printf("%f", n_dynamicmtx[i][j]);
    puts(seq1[0]);
    puts(seq2[0]);
    printf("alloclen = %d, headgp = %d, tailgp = %d\n", alloclen, headgp, tailgp);
    // Protein alphabets
    for(i = 0; i < nscoredalphabets; ++ i, putchar(' ')) printf("%c", amino[i]); puts("");
    // Gap char
	printf("Newgapstr = %s\n", newgapstr);
	//Gap penalty
	printf("Gap penalty = %d\n", ppenalty);
# endif

    /* Part 1: Defining the varibles of KBand */
    static TLS double **amino_dynamicmtx = NULL;
    static TLS int length1, length2, mxlength, band, needrerun;
	static TLS char *s1, *s2, *tmp, *res1, *res2, *tmp2;
	static TLS double oldval, val, gappenalty;
    
    /* Part 2: Allocing the varibles of KBand */
    amino_dynamicmtx = AllocateDoubleMtx(0x100, 0x100);

    /* Part 3: Give initial varibles of the varibles of Part 1 */
    //amino__dynamicmtx: weight of protein sequence
    for(i = 0; i < nalphabets; ++ i)
        for(j = 0; j < nalphabets; ++ j)
            amino_dynamicmtx[(unsigned char)amino[i]][(unsigned char)amino[j]] = (double)n_dynamicmtx[i][j];
    //lenght1, length2: length of sequence
    length1 = strlen(seq1[0]);
    length2 = strlen(seq2[0]);
	mxlength = MAX(length1, length2);
	oldval = -inf;
	s1 = seq1[0];
	s2 = seq2[0];
	band = alignband;
	gappenalty = (double)penalty;

    /* Part 4: Exception */
    // No data
    if(length1 == 0 || length2 == 0) 
	{
		FreeDoubleMtx(amino_dynamicmtx);
		return 0.0;
	}

    /* Part 5: KBand algorithm and alloc the matrix */
	res1 = AllocateCharVec((MAX(length1, length2) << 1) + 10);
	res2 = AllocateCharVec((MAX(length1, length2) << 1) + 10);


	if(alignband != NOTSPECIFIED) 
	{
		if(length1 < length2)
			val = Kband(s1, s2, length1, length2, gappenalty, band, amino_dynamicmtx, res1, res2, headgp, tailgp);
		else 
			val = Kband(s2, s1, length2, length1, gappenalty, band, amino_dynamicmtx, res2, res1, headgp, tailgp);
	}
	else
	{
	band = 10;
	if(length1 > length2)
	{
		// reporterr("Swapped\n");
		if(band > mxlength) val = Kband(s2, s1, length2, length1, gappenalty, band, amino_dynamicmtx, res2, res1, headgp, tailgp);
		while(band <= mxlength)
		{
			// According to the paper, band must be larger than abs(length1 - length2)
			val = Kband(s2, s1, length2, length1, gappenalty, band, amino_dynamicmtx, res2, res1, headgp, tailgp);
#if DEBUG
			reporterr("Score = %f, band = %d\n", val, band);
#endif
			if(val == -inf)
			{
				band <<= 1;
			}
			else if(val < oldval) {band >>= 1; needrerun = 1; break;}
			else if(val == oldval) {needrerun = 0; break;}
			else
			{
				oldval = val;
				band <<= 1;
				if(band > mxlength) break;
			}
		}
		// if(needrerun) Kband(s2, s1, length2, length1, gappenalty, band, amino_dynamicmtx, res2, res1, headgp, tailgp);
	}
	else
	{
		if(band > mxlength) val = Kband(s1, s2, length1, length2, gappenalty, band, amino_dynamicmtx, res1, res2, headgp, tailgp);
		while(band <= mxlength)
		{
			// According to the paper, band must be larger than abs(length1 - length2)
			val = Kband(s1, s2, length1, length2, gappenalty, band, amino_dynamicmtx, res1, res2, headgp, tailgp);
#if DEBUG
			reporterr("Score = %f, band = %d\n", val, band);
#endif
			if(val == -inf)
			{
				band <<= 1;
			}
			else if(val < oldval) {band >>= 1; needrerun = 1; break;}
			else if(val == oldval) {needrerun = 0; break;}
			else
			{
				oldval = val;
				band <<= 1;
				if(band > mxlength) break;
			}
		}
		// if(needrerun) Kband(s1, s2, length1, length2, gappenalty, band, amino_dynamicmtx, res1, res2, headgp, tailgp);
	}
	}

	
#if DEBUG
	printf("After alignment, s1 = %s, s2 = %s (Line 356 in Galign11.c)\n", res1, res2);
#endif

#if 1 // reverse tmp1, tmp2
	for(tmp2 = res1; *tmp2; ++ tmp2);
	for(tmp = s1; tmp2 != res1; ++ tmp) 
	{ 
		*tmp = *--tmp2;
		//printf("%s\n", s1);
	}
	*tmp = 0;
	for(tmp2 = res2; *tmp2; ++ tmp2);
	for(tmp = s2; tmp2 != res2; ++ tmp) 
	{
		//printf("%d %s\n", *tmp, res2);
		*tmp = *--tmp2;
		//printf("%s\n", s2);
	}
	*tmp = 0;
#else // copy directly
	strncpy(s1, res1, strlen(res1));
	strncpy(s2, res2, strlen(res2));
#endif

#ifndef NOTPRINTLEN
	if(strlen(s1) != strlen(res1))
	{
		reporterr("In G__align11, aligned length = (%d, %d)\n", strlen(s1), strlen(s2));
		if(strlen(s1) != strlen(s2)) reporterr("%s\n%s\n", s1, s2);
		reporterr("(%d, %d)\n", strlen(res1), strlen(res2));
	}
#endif

#if DEBUG
	printf("After process, s1 = %s, s2 = %s (Line 388 in Galign11.c)\n", s1, s2);
#endif

	/* Part 6: Free All varibles and return */
	FreeDoubleMtx(amino_dynamicmtx);
	free(res1);
	free(res2);
	return val;
}
