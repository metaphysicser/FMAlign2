#include "mltaln.h"
#include "dp.h"
#include "Kband.h"

#define DEBUG 0
#define protein_score(place1, place2) matrix_query(place1, place2, band, homo_matrix, len1 + 1, len2 + 1)

static const double inf = 1e100, eps = 1e-5;


void init__protein_score_Aalign(double **cpmx1, double **cpmx2, int len1, int len2, int band, double *ans)
{
	int i, j, k, l;
	double *scarr = AllocateDoubleVec(nalphabets + 5), *s, tmp;
	for(i = 0; i <= len1; ++ i)
	{
		for(s = scarr, j = 0; j < nalphabets; ++ j, ++ s) *s = cpmx1[j][i]; // NOTICE: cpmx[alphabet][place]
		for(j = left_band(len1 + 1, len2 + 1, i, band); j <= right_band(len1 + 1, len2 + 1, i, band); ++ j)	
		{
			for(s = scarr, k = 0; k < nalphabets; ++ k, ++ s)
				if (fabs(*s) > eps)
				{
					for(l = 0; l < nalphabets; ++ l)
					{
						if(fabs(tmp = cpmx2[l][j]) > eps)
						{
							// reporterr("(%d, %d) = (%c, %c), val = %f\n", i, j, amino[k], amino[l], amino_dis[(unsigned char)amino[k]][(unsigned char)amino[l]]);
							matrix_update(i, j, band, 
											*s * tmp * n_dis[k][l], 
											ans, len1 + 1, len2 + 1);
						}
					}
				}
		}
	}
	FreeDoubleVec(scarr);
	//fprintf(stderr, "OK protein score, band = %d\n", band);
}

double Kband__MSA_Aalign(int icyc, int jcyc, int len1, int len2, double **cpmx1, double **cpmx2, char *r1, char *r2, int band, double penalty)
{
	/*
		Every iteration of Kband.
		Return: Alignment score. -inf means this band can not finish this iteration.
	*/
	//fprintf(stderr, "In Kband, len1 = %d, len2 = %d, band = %d\n", len1, len2, band);
	//printf( "In Kband, len1 = %d, len2 = %d\n", len1, len2);
	int i, j, h, k;
	static TLS int len, swapped = 0;
	static TLS char *s;
	if(len1 > len2)
	{
		swapped = 1;
		len = len2; len2 = len1; len1 = len;
		s = r1; r1 = r2; r2 = s;
	}
	/* 
	   insert node: matrix_set(key, val)
	   find node: map_t *data = matrix_query(&tree, key)
	   The complexity of this block is : time O(n), space O(n)
	*/
	static TLS int max_place_i = 0, moveval, tmpval_col, tmpval_row;
	static TLS double score = 0, tmpval, max_data_i = -inf;
	static TLS char *tmpr1, *tmpr2;
	static TLS int *move_vec, *max_place_j, *mpjp;
	static TLS double *val_vec, *homo_matrix, *max_data_j, *mdjp;

	move_vec = AllocateIntVecLarge((long long)((band << 1 | 1) + (len2 - len1)) * (len2 + 10) + 10);
	val_vec = AllocateDoubleVecLarge((long long)((band << 1 | 1) + (len2 - len1)) * (len2 + 10) + 10);
	homo_matrix = AllocateDoubleVecLarge((long long)((band << 1 | 1) + (len2 - len1)) * (len2 + 10) + 10);
	max_place_j = AllocateIntVec(len2 + 10);
	max_data_j = AllocateDoubleVec(len2 + 10);


	//initial varibles
	init__protein_score_Aalign(cpmx1, cpmx2, len1, len2, band, homo_matrix);

	memset(r1, '?', sizeof(char) * (len1 + len2));
	memset(r2, '?', sizeof(char) * (len1 + len2));
	/* Initial of the matrix */
	matrix_set(0, 0, band, protein_score(0, 0), val_vec, len1 + 1, len2 + 1);
	for(i = 1; i <= len1; ++ i) matrix_set(i, 0, band, protein_score(i, 0), val_vec, len1 + 1, len2 + 1);
	for(i = 1; i <= len2; ++ i) matrix_set(0, i, band, protein_score(0, i), val_vec, len1 + 1, len2 + 1);
	/* Score Matrix and move_vec */
	max_data_i = -inf; max_place_i = 0; // Initial the max_data_i && max_place_i
	for(mdjp = max_data_j, mpjp = max_place_j, i = 0; i <= len2; ++ mdjp, ++ mpjp, ++ i) // Initial the max_data_j && max_place_j
	{
		*mdjp = matrix_query(0, i, band, val_vec, len1 + 1, len2 + 1);
		*mpjp = 0;
	}
	for(i = 1; i <= len1; ++ i)
	{
		max_data_i = matrix_query(i - 1, MAX(i - band, 0), band, val_vec, len1 + 1, len2 + 1) + penalty;
		max_place_i = MAX(i - band, 0);
		if(insideband(i - 1, 0, band, len1 + 1, len2 + 1))
		{
			tmpval = matrix_query(i - 1, 0, band, val_vec, len1 + 1, len2 + 1);
			if(*max_data_j < tmpval)
			{
				*max_data_j = tmpval;
				*max_place_j = i - 1;
			}
		}	
		for(j = left_band(len1 + 1, len2 + 1, i, band); j <= right_band(len1 + 1, len2 + 1, i, band); ++ j)
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
			if(i != len1 && j != len2) matrix_set(i, j, band, tmpval + protein_score(i, j), val_vec, len1 + 1, len2 + 1);
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
				if((tmpval_row = tmpval + penalty) >= max_data_i)
				{
					max_data_i = tmpval_row;
					max_place_i = j - 1;
					// reporterr("Max_data_i: (%d, %d) -> %d %f\n", i, j - 1, max_place_i, max_data_i);
				}
			}
			if(insideband(i - 1, j, band, len1 + 1, len2 + 1))
			{
				if((tmpval_col = tmpval + penalty) >= *(max_data_j + j))
				{
					*(max_data_j + j) = tmpval_col;
					*(max_place_j + j) = i - 1;
					// reporterr("Max_data_j: (%d, %d) -> %d %f\n", i - 1, j, *(max_place_j + j), *(max_data_j + j));
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

	/* Traceback */
	int limk = len1 + len2, nowi = len1, nowj = len2, val, tari, tarj;
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
		// reporterr("(nowi, nowj) = (%d, %d) -> (tari, tarj) = (%d, %d), k = %d, limk = %d, val = %d\n", nowi, nowj, tari, tarj, k, limk, val);
		// insert it!
		h = nowi - tari;
		while(-- h > 0) // if the difference is lower or same than 1, than DO NOT insert it!
		{
			*++ tmpr1 = 'o';
			*++ tmpr2 = *newgapstr;  
			++ k;
		}
		h = nowj - tarj;
		while(-- h > 0)
		{
			*++ tmpr1 = *newgapstr;
			*++ tmpr2 = 'o';
			++ k;
		}
		
		if(nowi <= 0 || nowj <= 0) break;
		*++ tmpr1 = 'o';
		*++ tmpr2 = 'o';
		++ k;
		nowi = tari, nowj = tarj;
		// reporterr("r1 = %s, r2 = %s\n", r1, r2);
	}
	*++ tmpr1 = 0;
	*++ tmpr2 = 0;

	score = matrix_query(len1, len2, band, val_vec, len1 + 1, len2 + 1);
	FreeDoubleVec(val_vec);
	FreeDoubleVec(homo_matrix);
	FreeIntVec(move_vec);
	FreeIntVec(max_place_j);
	FreeDoubleVec(max_data_j);
#if DEBUG
	printf("%s %s\n", r1, r2);
#endif
	if(swapped)
	{
		len = len2; len2 = len1; len1 = len;
		s = r1; r1 = r2; r2 = s;
	}

	return score;

}

double Aalign( char **seq1, char **seq2, double *eff1, double *eff2, int icyc, int jcyc, int alloclen )
{
    /* 
        Kband Algorithm: align **two** **Profile** sequences
        Return: Alignment Score. The sequence seq1 & seq2 must be aligned.
    */
    // Part -1: loop varibles
    register int i, j;
    /* Part 1: Defining the varibles of KBand */
    static TLS double **amino_dynamicmtx = NULL;
    static TLS int length1, length2, tmplen, band, swapped = 0, swapp, needrerun = 0, mxlength;
	static TLS char **res1, **res2; // Result array
	static TLS char *gap1, *gap2, *swapgap; // Gap string '-' has gap 'o' not gap
	static TLS int gaplen1, gaplen2; // Gap length without '?'
	static TLS double oldval, val, gappenalty;
	static TLS double **cpmx1, **cpmx2; // Every amino in every place
	char *rev, *rev2;

    /* Part 3: Give initial varibles of the varibles of Part 1 */
	band = alignband;
	gappenalty = (double)penalty * 0.5;
    length1 = strlen(seq1[0]);
    length2 = strlen(seq2[0]);
	mxlength = MAX(length1, length2);
	gaplen1 = 0;
	gaplen2 = 0;
    /* Part 4: Exception */
    // No data
    if(length1 == 0 || length2 == 0) return 0.0;
	
	/* Part 2: Allocing the varibles that using length information */
	length1 += 10;
	length2 += 10;

	cpmx1 = AllocateDoubleMtx(nalphabets, length1);
	cpmx2 = AllocateDoubleMtx(nalphabets, length2);
	gap1 = AllocateCharVec(length1 + length2);
	gap2 = AllocateCharVec(length1 + length2);
	swapgap = AllocateCharVec(length1 + length2);
	res1 = AllocateCharMtx(icyc, length1 + length2 + 80);
	res2 = AllocateCharMtx(jcyc, length1 + length2 + 80);

	length1 -= 10;
	length2 -= 10;

    /* Part 5: KBand algorithm and alloc the matrix */
	cpmx_calc_new(seq1, cpmx1, eff1, length1, icyc);
	cpmx_calc_new(seq2, cpmx2, eff2, length2, jcyc);
	if(alignband != NOTSPECIFIED)
	{
		if(length1 < length2) 
			val = Kband__MSA_Aalign(icyc, jcyc, length1, length2, cpmx1, cpmx2, gap1, gap2, band, gappenalty);
		else
			val = Kband__MSA_Aalign(jcyc, icyc, length2, length1, cpmx2, cpmx1, gap2, gap1, band, gappenalty);
	}
	else
	{
	band = 10;
	if(length1 < length2)
	{
		if(band > mxlength) val = Kband__MSA_Aalign(icyc, jcyc, length1, length2, cpmx1, cpmx2, gap1, gap2, band, gappenalty);
		while(band <= mxlength)
		{
			val = Kband__MSA_Aalign(icyc, jcyc, length1, length2, cpmx1, cpmx2, gap1, gap2, band, gappenalty);
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
		// if(needrerun) Kband__MSA_Aalign(icyc, jcyc, length1, length2, cpmx1, cpmx2, gap1, gap2, band, gappenalty);
	}
	else
	{
		if(band > mxlength) val = Kband__MSA_Aalign(jcyc, icyc, length2, length1, cpmx2, cpmx1, gap2, gap1, band, gappenalty);
		while(band <= mxlength)
		{
			val = Kband__MSA_Aalign(jcyc, icyc, length2, length1, cpmx2, cpmx1, gap2, gap1, band, gappenalty);
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
		// if(needrerun) Kband__MSA_Aalign(jcyc, icyc, length2, length1, cpmx2, cpmx1, gap2, gap1, band, gappenalty);
	}
	}

	for(rev2 = gap1; ((*rev2) && *rev2 != '?'); ++ rev2) ++ gaplen1;
	for(rev = swapgap; rev2 != gap1; ++ rev) 
	{ 
		*rev = *--rev2;
		//reporterr("%s\n", swapgap);
	}
	*rev = '\0';
	strncpy(gap1, swapgap, gaplen1 * sizeof(char));
	for(rev2 = gap2; ((*rev2) && *rev2 != '?'); ++ rev2) ++ gaplen2;
	for(rev = swapgap; rev2 != gap2; ++ rev) 
	{
		//printf("%d %s\n", *tmp, res2);
		*rev = *--rev2;
		//reporterr("%s\n", swapgap);
	}
	*rev = '\0';
	strncpy(gap2, swapgap, gaplen2 * sizeof(char));

	/* According to gap1 & gap2, insert the gap into the profiles */
	if(strchr(gap1, '-'))
	{
		for(i = 0; i < icyc; ++ i) gapireru(res1[i], seq1[i], gap1);
		for(i = 0; i < icyc; ++ i) strcpy(seq1[i], res1[i]);
	}
	if(strchr(gap2, '-'))
	{
		for(j = 0; j < jcyc; ++ j) gapireru(res2[j], seq2[j], gap2);
		for(j = 0; j < jcyc; ++ j) strcpy(seq2[j], res2[j]);
	}

	FreeDoubleMtx(cpmx1);
	FreeDoubleMtx(cpmx2);
	free(gap1);
	free(gap2);
	free(swapgap);
	FreeCharMtx(res1);
	FreeCharMtx(res2);
	return val;
}
