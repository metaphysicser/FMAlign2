#include "mltaln.h"
#include "dp.h"
#include "Kband.h"

#define USEDEXGAP 1
#define NOTPRINT
#define DEBUG 0
#define NOTPRINTLEN

#define protein_score(place1, place2) matrix_query(place1, place2, band, homo_matrix, len1 + 1, len2 + 1)

static const double inf = 1e100, eps = 1e-5;

void init__protein_score(double **cpmx1, double **cpmx2, int len1, int len2, int band, double **amino_dis, double *ans)
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
											*s * tmp * amino_dis[(unsigned char)amino[k]][(unsigned char)amino[l]], 
											ans, len1 + 1, len2 + 1);
						}
					}
				}
		}
	}
	FreeDoubleVec(scarr);
	//fprintf(stderr, "OK protein score, band = %d\n", band);
#if DEBUG
	for(i = 0; i <= len1; ++ i)
		for(j = 0; j <= len2; ++ j)
		{
			if(insideband(i, j, band, len1 + 1, len2 + 1)) 
			{
				reporterr("(%d %d): %f\n", i, j, matrix_query(i, j, band, ans, len1 + 1, len2 + 1));
				for(k = 0; k < nalphabets; ++ k) 
				{
					if(fabs(cpmx1[k][i]) > eps) 
					{
						reporterr("cpmx1[%d][%d] = %f\n", k, i, cpmx1[k][i]);
						for(l = 0; l < nalphabets; ++ l) 
							if(fabs(cpmx2[l][j]) > eps) 
								reporterr("cpmx2[%d][%d] = %f, number = %f\n", l, j, cpmx2[l][j], amino_dis[(unsigned char)amino[k]][(unsigned char)amino[l]]);
					}
				}
				 
			}
		}
#endif
}


double Kband__MSA(int p1, int p2, int len1, int len2, int band, 
                  double *og1, double *fg1, double *og2, double *fg2, double *gf1, double *gf2, double **cpmx1, double **cpmx2,
				  char *r1, char *r2, double **ad, double hgp1, double hgp2, int headgp, int tailgp)
{
	/*
		Every iteration of Kband.
		Return: Alignment score. -inf means this band can not finish this iteration.
	*/
# if DEBUG
	fprintf(stderr, "In Kband, len1 = %d, len2 = %d, band = %d\n", len1, len2, band);
	printf( "In Kband, len1 = %d, len2 = %d\n", len1, len2);
# endif
	int i, j, k, len = MAX(len1, len2);
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
	init__protein_score(cpmx1, cpmx2, len1, len2, band, ad, homo_matrix);

	memset(r1, '?', sizeof(char) * (len1 + len2));
	memset(r2, '?', sizeof(char) * (len1 + len2));
	/* Initial of the matrix */
	matrix_set(0, 0, band, protein_score(0, 0), val_vec, len1 + 1, len2 + 1);
	for(i = 1; i <= len1; ++ i) matrix_set(i, 0, band, protein_score(i, 0), val_vec, len1 + 1, len2 + 1);
	for(i = 1; i <= len2; ++ i) matrix_set(0, i, band, protein_score(0, i), val_vec, len1 + 1, len2 + 1);
	/* headgap may be 0? */
	if(headgp == 1)
	{
		for(i = 1; i <= len1; ++ i) matrix_update(i, 0, band, *og1 * hgp2 + *(fg1 + i - 1) * *gf2, val_vec, len1 + 1, len2 + 1);
		for(i = 1; i <= len2; ++ i) matrix_update(0, i, band, *og2 * hgp1 + *(fg2 + i - 1) * *gf1, val_vec, len1 + 1, len2 + 1);
	}
	
	/* Score Matrix and move_vec */
	max_data_i = -inf; max_place_i = 0; // Initial the max_data_i && max_place_i
	for(mdjp = max_data_j, mpjp = max_place_j, i = 0; i <= len2; ++ mdjp, ++ mpjp, ++ i) // Initial the max_data_j && max_place_j
	{
		*mdjp = matrix_query(0, i, band, val_vec, len1 + 1, len2 + 1) + *(og1 + 1) * *(gf2 + i - 1);
		*mpjp = 0;
	}
	for(i = 1; i <= len1; ++ i)
	{
		max_data_i = matrix_query(i - 1, left_band(len1 + 1, len2 + 1, i - 1, band), band, val_vec, len1 + 1, len2 + 1) + *(og2 + left_band(len1 + 1, len2 + 1, i - 1, band) + 1) * *(gf1 + i - 1);
		max_place_i = left_band(len1 + 1, len2 + 1, i - 1, band);
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
				if((tmpval_row = max_data_i + *(fg2 + j - 1) * *(gf1 + i)) > tmpval)
				{
					tmpval = tmpval_row;
					moveval = -(j - max_place_i);
					// reporterr("Set by i, j - 1, j = %d, max_place_i = %d, val = %f\n", j, max_place_i, tmpval);
				}
				
			}
			if(insideband(i - 1, j, band, len1 + 1, len2 + 1))
			{
				if((tmpval_col = *(max_data_j + j) + *(gf2 + j) * *(fg1 + i - 1)) > tmpval)
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
				if((tmpval_row = tmpval + *(og2 + j) * *(gf1 + i - 1)) >= max_data_i)
				{
					max_data_i = tmpval_row;
					max_place_i = j - 1;
					// reporterr("Max_data_i: (%d, %d) -> %d %f\n", i, j - 1, max_place_i, max_data_i);
				}
			}
			if(insideband(i - 1, j, band, len1 + 1, len2 + 1))
			{
				if((tmpval_col = tmpval + *(og1 + i) * *(gf2 + j - 1)) >= *(max_data_j + j))
				{
					*(max_data_j + j) = tmpval_col;
					*(max_place_j + j) = i - 1;
					// reporterr("Max_data_j: (%d, %d) -> %d %f\n", i - 1, j, *(max_place_j + j), *(max_data_j + j));
				}
			}
		}

		/* assure that the place appears in max_place_j is in the band */
		for(k = 0, mdjp = max_data_j, mpjp = max_place_j; k < left_band(len1 + 1, len2 + 1, i - 1, band); ++ k, ++ mdjp, ++ mpjp) // clean the outside band date (only need to it?)
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
		j = nowi - tari;
		while(-- j > 0) // if the difference is lower or same than 1, than DO NOT insert it!
		{
			*++ tmpr1 = 'o';
			*++ tmpr2 = *newgapstr;  
			++ k;
		}
		j = nowj - tarj;
		while(-- j > 0)
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
	return score;
}

double A__align( double **n_dynamicmtx, int penalty_l, int penalty_ex_l, char **seq1, char **seq2, double *eff1, double *eff2, int icyc, int jcyc, int alloclen, char *sgap1, char *sgap2, char *egap1, char *egap2, int headgp, int tailgp)
{
    //puts("In A__align11");
    /* 
        Kband Algorithm: align **two** **Profile** sequences
        Return: Alignment Score. The sequence seq1 & seq2 must be aligned.
    */
    // Part -1: loop varibles
    int i, j;
# ifndef NOTPRINT
    // Part 0: Print some arguments
    printf("nalphabets = %d, nscoredalphabets = %d\n", nalphabets, nscoredalphabets);
    for(i = 0; i < nscoredalphabets; ++ i, puts(""))
        for(j = 0; j < nscoredalphabets; ++ j, putchar(' '))
            printf("%f", n_dynamicmtx[i][j]);
    for(i = 0; i < icyc; ++ i) puts(seq1[i]);
	puts("");
    for(j = 0; j < jcyc; ++ j) puts(seq2[j]);
	// the number of Sequences
	printf("icyc = %d, jcyc = %d\n", icyc, jcyc);
    printf("alloclen = %d, headgp = %d, tailgp = %d\n", alloclen, headgp, tailgp);
    // Protein alphabets
    for(i = 0; i < nscoredalphabets; ++ i, putchar(' ')) printf("%c", amino[i]); puts("");
    // Gap char
	printf("Newgapstr = %s\n", newgapstr);
	// Gap penalty
	printf("Gap penalty = %d, penalty_ex = %d\n", penalty_l, penalty_ex_l);
	// Start gap
	printf("Start gap pointer = %p, %p\n", sgap1, sgap2);
	// End gap
	printf("End gap pointer = %p, %p\n", egap1, egap2);
# endif

    /* Part 1: Defining the varibles of KBand */
    static TLS double **amino_dynamicmtx = NULL;
    static TLS int length1, length2, mxlength, tmplen, band, swapped = 0, swapp, needrerun;
	static TLS char **res1, **res2; // Result array
	static TLS char *gap1, *gap2, *swapgap; // Gap string '-' has gap 'o' not gap
	static TLS int gaplen1, gaplen2; // Gap length without '?'
	static TLS double oldval, val, gappenalty;
	static TLS double **cpmx1, **cpmx2; // Every amino in every place
	static TLS double *gapfreq1, *gapfreq2; // Gap frequency
	double *gapf1qp, *gapf2qp;
	char *rev, *rev2;
#if USEDEXGAP
	static TLS double *opgap1f, *opgap2f; // Opening Gap frequency
	static TLS double *fgap1f, *fgap2f; // Closing gap frequency
    double *opgap1fp, *opgap2fp, *fgap1fp, *fgap2fp;
	static TLS double headgapfreq1, headgapfreq2; // first place penalty
#endif

    /* Part 2: Allocing the varible of KBand - amino matrix */
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
	band = alignband;
	gappenalty = (double)penalty;

    /* Part 4: Exception */
    // No data
    if(length1 == 0 || length2 == 0) 
	{
		FreeDoubleMtx(amino_dynamicmtx);
		return 0.0;
	}

	/* Part 2.2: Allocing the varibles that using length information */
	length1 += 10;
	length2 += 10;

	cpmx1 = AllocateDoubleMtx(nalphabets, length1);
	cpmx2 = AllocateDoubleMtx(nalphabets, length2);
	gapfreq1 = AllocateDoubleVec(length1);gapf1qp = gapfreq1;
	gapfreq2 = AllocateDoubleVec(length2);gapf2qp = gapfreq2;

	opgap1f = AllocateDoubleVec(length1); opgap1fp = opgap1f;
	opgap2f = AllocateDoubleVec(length2); opgap2fp = opgap2f;
	fgap1f = AllocateDoubleVec(length1);  fgap1fp = fgap1f;
	fgap2f = AllocateDoubleVec(length2);  fgap2fp = fgap2f;

	res1 = AllocateCharMtx(icyc, length1 + length2 + 80);
	res2 = AllocateCharMtx(jcyc, length1 + length2 + 80);
	gap1 = AllocateCharVec(length1 + length2);
	gap2 = AllocateCharVec(length1 + length2);
	swapgap = AllocateCharVec(length1 + length2);

	length1 -= 10;
	length2 -= 10;

	gaplen1 = 0;
	gaplen2 = 0;
    /* Part 5: KBand algorithm and alloc the matrix */

	// Calcuate the frequency of all characters
	cpmx_calc_new(seq1, cpmx1, eff1, length1, icyc);
#ifndef NOTPRINTLEN
	fprintf(stderr, "OK on cpmx_calc to seq1\n");
#endif
	cpmx_calc_new(seq2, cpmx2, eff2, length2, jcyc);
#ifndef NOTPRINTLEN
	fprintf(stderr, "OK on cpmx_calc to seq2\n");
#endif
	// Count open gap && final gap
	st_OpeningGapCount(opgap1fp, icyc, seq1, eff1, length1);
	st_FinalGapCount(fgap1fp, icyc, seq1, eff1, length1);

	st_OpeningGapCount(opgap2fp, jcyc, seq2, eff2, length2);
	st_FinalGapCount(fgap2fp, jcyc, seq2, eff2, length2);
	
	// get head gap frequency && tail gap frequency
	outgapcount(&headgapfreq1, icyc, sgap1, eff1);
	outgapcount(&headgapfreq2, jcyc, sgap2, eff2);

	outgapcount( gapf1qp + length1, icyc, egap1, eff1 );
	outgapcount( gapf2qp + length2, jcyc, egap2, eff2 );

	if( legacygapcost == 0 )
	{
		// Count common gap frequency 
		gapcountf(gapf1qp, seq1, icyc, eff1, length1);
		gapcountf(gapf2qp, seq2, jcyc, eff2, length2);

		// Reverse
		for(i = 0; i <= length1; ++ i) gapf1qp[i] = 1.0 - gapf1qp[i];
		for(i = 0; i <= length2; ++ i) gapf2qp[i] = 1.0 - gapf2qp[i];
		headgapfreq1 = 1.0 - headgapfreq1;
		headgapfreq2 = 1.0 - headgapfreq2;
	}
	else
	{
		for(i = 0; i <= length1; ++ i) gapf1qp[i] = 1.0;
		for(i = 0; i <= length2; ++ i) gapf2qp[i] = 1.0;
		headgapfreq1 = 1.0;
		headgapfreq2 = 1.0;
	}

	// Calcuate the opening gap && final gap
	for(i = 0; i <= length1; ++ i ) 
	{
		opgap1f[i] = 0.5 * ( 1.0 - opgap1f[i] ) * gappenalty * ( gapf1qp[i] );
		fgap1f[i] = 0.5 * ( 1.0 - fgap1f[i] ) * gappenalty * ( gapf1qp[i] );
	}
	for(i = 0; i <= length2; ++ i ) 
	{
		opgap2f[i] = 0.5 * ( 1.0 - opgap2f[i] ) * gappenalty * ( gapf2qp[i] );
		fgap2f[i] = 0.5 * ( 1.0 - fgap2f[i] ) * gappenalty * ( gapf2qp[i] );
	}


#if 0
	reporterr("\nopening gap1: \n");
	for(i = 0; i <= length1; ++ i, reporterr(" ")) reporterr("%f", *(opgap1f + i));
	reporterr("\nopening gap2: \n");
	for(i = 0; i <= length2; ++ i, reporterr(" ")) reporterr("%f", *(opgap2f + i));
	reporterr("\nfinal gap1: \n");
	for(i = 0; i <= length1; ++ i, reporterr(" ")) reporterr("%f", *(fgap1f + i));
	reporterr("\nfinal gap2: \n");
	for(i = 0; i <= length2; ++ i, reporterr(" ")) reporterr("%f", *(fgap2f + i));
	reporterr("\ngap1 penalty: \n");
	for(i = 0; i <= length1; ++ i, reporterr(" ")) reporterr("%f", *(gapf1qp + i));
	reporterr("\ngap2 penalty: \n");
	for(i = 0; i <= length2; ++ i, reporterr(" ")) reporterr("%f", *(gapf2qp + i));
	reporterr("\n");
#endif
	if(alignband != NOTSPECIFIED)
	{
		if(length1 < length2) 
			val = Kband__MSA(icyc, jcyc, length1, length2, band, opgap1f, fgap1f, opgap2f, fgap2f, gapf1qp, gapf2qp, 
		                     cpmx1, cpmx2, gap1, gap2, amino_dynamicmtx, headgapfreq1, headgapfreq2, headgp, tailgp);
		else
			val = Kband__MSA(jcyc, icyc, length2, length1, band, opgap2f, fgap2f, opgap1f, fgap1f, gapf2qp, gapf1qp, 
			                 cpmx2, cpmx1, gap2, gap1, amino_dynamicmtx, headgapfreq2, headgapfreq1, headgp, tailgp);
		
	}
	else
	{
	band = 10;
	if(length1 < length2)
	{
		val = Kband__MSA(icyc, jcyc, length1, length2, band, opgap1f, fgap1f, opgap2f, fgap2f, gapf1qp, gapf2qp, 
			             cpmx1, cpmx2, gap1, gap2, amino_dynamicmtx, headgapfreq1, headgapfreq2, headgp, tailgp);
		while(band <= mxlength)
		{
		// According to the paper, band must be larger than abs(length1 - length2)
			val = Kband__MSA(icyc, jcyc, length1, length2, band, opgap1f, fgap1f, opgap2f, fgap2f, gapf1qp, gapf2qp, 
			                 cpmx1, cpmx2, gap1, gap2, amino_dynamicmtx, headgapfreq1, headgapfreq2, headgp, tailgp);
#if DEBUG
			reporterr("Val = %f, band = %d, length1 = %d, length2 = %d\n", val, band, length1, length2);
#endif
			if(val == -inf)
			{
				band <<= 1;
			}
			else if(val == oldval) {band >>= 1; needrerun = 1; break;}
			else
			{
				oldval = val;
				band <<= 1;
				if(band > mxlength) break;
			}
		}
	}
	else
	{
		val = Kband__MSA(jcyc, icyc, length2, length1, band, opgap2f, fgap2f, opgap1f, fgap1f, gapf2qp, gapf1qp, 
			             cpmx2, cpmx1, gap2, gap1, amino_dynamicmtx, headgapfreq2, headgapfreq1, headgp, tailgp);
		while(band <= mxlength)
		{
		// According to the paper, band must be larger than abs(length1 - length2)
			val = Kband__MSA(jcyc, icyc, length2, length1, band, opgap2f, fgap2f, opgap1f, fgap1f, gapf2qp, gapf1qp, 
			                 cpmx2, cpmx1, gap2, gap1, amino_dynamicmtx, headgapfreq2, headgapfreq1, headgp, tailgp);
#if DEBUG
			reporterr("Val = %f, band = %d, length1 = %d, length2 = %d\n", val, band, length1, length2);
#endif
			if(val == -inf)
			{
				band <<= 1;
			}
			else if(val == oldval) {band >>= 1; needrerun = 1; break;}
			else
			{
				oldval = val;
				band <<= 1;
				if(band > mxlength) break;
			}
		}
	}
	}

#if DEBUG
	reporterr("%s\n", gap1);
#endif
	for(rev2 = gap1; ((*rev2) && *rev2 != '?'); ++ rev2) ++ gaplen1;
	for(rev = swapgap; rev2 != gap1; ++ rev) 
	{ 
		*rev = *--rev2;
		//reporterr("%s\n", swapgap);
	}
	*rev = '\0';
#if DEBUG
	reporterr("gaplen = %d\n", gaplen1);
#endif
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


	/* Part 6: Free All varibles and return */
	FreeDoubleMtx(amino_dynamicmtx);
	FreeCharMtx(res1);
	FreeCharMtx(res2);
	FreeDoubleMtx(cpmx1);
	FreeDoubleMtx(cpmx2);
	free(swapgap);
	free(gap1);
	free(gap2);
	FreeDoubleVec(gapfreq1);
	FreeDoubleVec(gapfreq2);
#if USEDEXGAP
	FreeDoubleVec(fgap1f);
	FreeDoubleVec(fgap2f);
	FreeDoubleVec(opgap1f);
	FreeDoubleVec(opgap2f);
#endif
#ifndef NOTPRINTLEN
	reporterr("In A__align, [");
	for(i = 0; i < icyc; ++ i) reporterr("%d ", strlen(seq1[i]));
	reporterr(",");
	for(j = 0; j < jcyc; ++ j) reporterr("%d ", strlen(seq2[j]));
	reporterr("]\n");
#endif
	return val;
}
