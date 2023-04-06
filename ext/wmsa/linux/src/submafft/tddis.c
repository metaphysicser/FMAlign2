#include "mltaln.h"

#define DEBUG 0
		
void cpmx_calc( char **seq, double **cpmx, double *eff, int lgth, int clus )
{
	int  i, j, k;
	double totaleff = 0.0;

	for( i=0; i<clus; i++ ) totaleff += eff[i]; 
	for( i=0; i<nalphabets; i++ ) for( j=0; j<lgth; j++ ) cpmx[i][j] = 0.0;
	for( j=0; j<lgth; j++ ) for( k=0; k<clus; k++ )
			cpmx[(int)amino_n[(unsigned char)seq[k][j]]][j] += (double)eff[k] / totaleff;
}


int fastconjuction_noname( int *memlist, char **seq, char **aseq, double *peff, double *eff, char *d, double mineff, double *oritotal  )
{
	int m, k, dln;
	char b[B];
	double total;

#if DEBUG
	fprintf( stderr, "s = %d\n", s );
#endif

	total = 0.0;
	d[0] = 0;
	dln = 0;
	for( k=0; *memlist!=-1; memlist++, k++ )
	{
		m = *memlist;
		dln += sprintf( b, " %d", m+1 ); 
#if 1
		if( dln < 100 ) strcat( d, b );
#else
		strcat( d, b );
#endif
		aseq[k] = seq[m];
		if( eff[m] < mineff )
			peff[k] = mineff;
		else
			peff[k] = eff[m];

		total += peff[k];
	}
	if( oritotal ) *oritotal = total;
#if 1
	for( m=0; m<k; m++ )
	{
//		fprintf( stderr, "Apr17   peff[%d] = %20.10f\n", m, peff[m] );
		peff[m] /= total;
	}
#endif
	return( k );
}
void cpmx_calc_new( char **seq, double **cpmx, double *eff, int lgth, int clus ) // summ eff must be 1.0
{
	int  i, j, k;
	double feff;
	double *cpmxpt, **cpmxptpt;
	char *seqpt;

	j = nalphabets;
	cpmxptpt = cpmx;
	while( j-- )
	{
		cpmxpt = *cpmxptpt++;
		i = lgth;
		while( i-- )
			*cpmxpt++ = 0.0;
	}
	for( k=0; k<clus; k++ )
	{
		feff = (double)eff[k];
		seqpt = seq[k];
		// fprintf( stderr, "seqpt = %s, seqlen = %d, lgth = %d\n", seqpt, strlen(seqpt), lgth );
		for( j=0; j<lgth; j++ )
		{
			cpmx[(unsigned char)amino_n[(unsigned char)*seqpt++]][j] += feff;
		}
	}
}
void MScpmx_calc_new( char **seq, double **cpmx, double *eff, int lgth, int clus ) // summ eff must be 1.0
{
	int  i, j, k;
	double feff;
	double **cpmxptpt, *cpmxpt;
	char *seqpt;

	j = lgth;
	cpmxptpt = cpmx;
	while( j-- )
	{
		cpmxpt = *cpmxptpt++;
		i = nalphabets;
		while( i-- )
			*cpmxpt++ = 0.0;
	}
	for( k=0; k<clus; k++ )
	{
		feff = (double)eff[k];
		seqpt = seq[k];
		cpmxptpt = cpmx;
		j = lgth;
		while( j-- )
			(*cpmxptpt++)[(int)amino_n[(unsigned char)*seqpt++]] += feff;
	}
#if 0
	for( j=0; j<lgth; j++ ) for( i=0; i<nalphabets; i++ ) cpmx[j][i] = 0.0;
	for( k=0; k<clus; k++ )
	{
		feff = (double)eff[k];
		for( j=0; j<lgth; j++ ) 
			cpmx[j][(int)amino_n[(int)seq[k][j]]] += feff;
	}
#endif
}

void cpmx_calc_add( char **seq, double **cpmx, double *eff, int lgth, int clus ) // lastmem = newmem; summ eff must be 1.0
{
	double neweff, orieff;
	int newmem, i, j;

	newmem = clus-1;
	neweff = eff[clus-1];
	orieff = 1.0-neweff;
#if 1 // TESTING Feb/1/18:00
	for( j=0; j<lgth; j++ )
	{
		for( i=0;i<nalphabets; i++ ) cpmx[i][j] *= orieff;
		cpmx[(unsigned char)amino_n[(unsigned char)seq[newmem][j]]][j] += neweff;
	}
#else // possibly faster?
	for( i=0;i<nalphabets; i++ ) 
	{
		for( j=0; j<lgth; j++ ) cpmx[i][j] *= orieff;
	}
	for( j=0; j<lgth; j++ ) cpmx[(unsigned char)amino_n[(unsigned char)seq[newmem][j]]][j] += neweff;
#endif
}
