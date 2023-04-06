#include "mltaln.h"
static TLS int n20or4or2;
#define DEBUG 0

#if DEBUG
static FILE *fftfp;
#endif

static void vec_init( Fukusosuu *result, int nlen )
{
	while( nlen-- )
	{
		result->R = result->I = 0.0;
		result++;
	}
}

static void seq_vec_2( Fukusosuu *result, double *score, double incr, char *seq )
{
	static TLS int n;
	for( ; *seq; result++ )
	{
		n = amino_n[(int)*seq++];
		if( n < 20 && n >= 0 ) result->R += incr * score[n];
#if 0
		fprintf( stderr, "n=%d, score=%f, inc=%f R=%f\n",n,  score[n], incr * score[n], result->R );
#endif
	}
}

static void seq_vec_4( Fukusosuu *result, double incr, char *seq )
{
	char s;
	for( ; *seq; result++ )
	{
		s = *seq++;
		if( s == 'a' )
			result->R += incr;
		else if( s == 't' )
			result->R -= incr;
		else if( s == 'g' )
			result->I += incr;
		else if( s == 'c' )
			result->I -= incr;
	}
}

static void seq_vec_3( Fukusosuu **result, double incr, char *seq )
{
	int i;
	int n;
	for( i=0; *seq; i++ )
	{
		n = amino_n[(int)*seq++];
		if( n < n20or4or2 && n >= 0 ) result[n][i].R += incr;
	}
}

static void seq_vec_5( Fukusosuu *result, double *score1, double *score2, double incr, char *seq )
{
	int n;
	for( ; *seq; result++ )
	{
		n = amino_n[(int)*seq++];
		if( n > 20 ) continue;
		result->R += incr * score1[n];
		result->I += incr * score2[n];
#if 0
		fprintf( stderr, "n=%d, score=%f, inc=%f R=%f\n",n,  score[n], incr * score[n], result->R );
#endif
	}
}

static void mymergesort( int first, int last, Segment **seg )
{
	int middle;
	static TLS int i, j, k, p;
	static TLS int allo = 0;
	static TLS Segment **work = NULL;

	if( seg == NULL )
	{
		if( work ) free( work ); 
		work = NULL;
		allo = 0;
		return;
	}

	if( last > allo )
	{
		allo = last;
		if( work ) free( work );
		work = (Segment **)calloc( allo / 2 + 1, sizeof( Segment *) );
	}

	if( first < last )
	{
		middle = ( first + last ) / 2;
		mymergesort( first, middle, seg );
		mymergesort( middle+1, last, seg );
		p = 0;
		for( i=first; i<=middle; i++ ) work[p++] = seg[i];
		i = middle + 1; j = 0; k = first;
		while( i <= last && j < p )
		{
			if( work[j]->center <= seg[i]->center ) 
				seg[k++] = work[j++];
			else
				seg[k++] = seg[i++];
		}
		while( j < p ) seg[k++] = work[j++];
	}
}

double Falign( int **whichmtx, double ***scoringmatrices, double **n_dynamicmtx,
			  char  **seq1, char  **seq2, 
			  double *eff1, double *eff2, 
			  double **eff1s, double **eff2s,
			  int    clus1, int    clus2,
			  int alloclen, int *fftlog,
			  int *chudanpt, int chudanref, int *chudanres )
{
	int i, j, k, l, m, maxk;
	int nlen, nlen2, nlen4;
	static TLS int crossscoresize = 0;
	char **tmpseq1 = NULL;
	char **tmpseq2 = NULL;
	char **tmpptr1 = NULL;
	char **tmpptr2 = NULL;
	char **tmpres1 = NULL;
	char **tmpres2 = NULL;
	char **result1 = NULL;
	char **result2 = NULL;

	static TLS Fukusosuu **seqVector1 = NULL;
	static TLS Fukusosuu **seqVector2 = NULL;
	static TLS Fukusosuu **naiseki = NULL;   
	static TLS Fukusosuu *naisekiNoWa = NULL; 
	static TLS double *soukan = NULL;
	static TLS double **crossscore = NULL;
	int nlentmp;
	static TLS int *kouho = NULL;
	static TLS Segment *segment = NULL;
	static TLS Segment *segment1 = NULL;
	static TLS Segment *segment2 = NULL;
	static TLS Segment **sortedseg1 = NULL;
	static TLS Segment **sortedseg2 = NULL;
	static TLS int *cut1 = NULL;
	static TLS int *cut2 = NULL;
	char *sgap1, *egap1, *sgap2, *egap2;
	static TLS int localalloclen = 0;
	int lag;
	int tmpint;
	int count, count0;
	int len1, len2;
	int totallen;
	double totalscore;
	double dumdb = 0.0;
	int headgp, tailgp;

	len1 = strlen( seq1[0] );
	len2 = strlen( seq2[0] );
	nlentmp = MAX( len1, len2 );

	nlen = 1;
	while( nlentmp >= nlen ) nlen <<= 1;

	nlen2 = nlen/2; nlen4 = nlen2 / 2;


	result1 = AllocateCharMtx( clus1, alloclen );
	result2 = AllocateCharMtx( clus2, alloclen );
	tmpres1 = AllocateCharMtx( clus1, alloclen ); //reporterr("%d %d %p\n", clus1, alloclen, tmpres1);
	tmpres2 = AllocateCharMtx( clus2, alloclen ); //reporterr("%d %d %p\n", clus2, alloclen, tmpres2);
	sgap1 = AllocateCharVec( clus1 );
	egap1 = AllocateCharVec( clus1 );
	sgap2 = AllocateCharVec( clus2 );
	egap2 = AllocateCharVec( clus2 );
	tmpptr1 = AllocateCharMtx( clus1, 0 );
	tmpptr2 = AllocateCharMtx( clus2, 0 );
	tmpseq1 = AllocateCharMtx( clus1, nlen );
	tmpseq2 = AllocateCharMtx( clus2, nlen );
	//++ tmpres1; ++ tmpres2;
	if( !localalloclen )
	{
		kouho = AllocateIntVec( NKOUHO );
		cut1 = AllocateIntVec( MAXSEG );
		cut2 = AllocateIntVec( MAXSEG );
//		crossscore = AllocateDoubleMtx( MAXSEG, MAXSEG );
		segment = (Segment *)calloc( MAXSEG, sizeof( Segment ) );
		segment1 = (Segment *)calloc( MAXSEG, sizeof( Segment ) );
		segment2 = (Segment *)calloc( MAXSEG, sizeof( Segment ) );
		sortedseg1 = (Segment **)calloc( MAXSEG, sizeof( Segment * ) );
		sortedseg2 = (Segment **)calloc( MAXSEG, sizeof( Segment * ) );
		if( !( segment && segment1 && segment2 && sortedseg1 && sortedseg2 ) )
			ErrorExit( "Allocation error\n" );

		if     ( scoremtx == -1 ) n20or4or2 = 1;
		else if( fftscore )       n20or4or2 = 1;
		else                      n20or4or2 = 20;
	}
#if DEBUG
	reporterr("len1 = %d, len2 = %d\n", len1, len2 );
	reporterr("nlentmp = %d, nlen = %d\nscoremtx = %d, fftscore = %d\n", nlentmp, nlen, scoremtx, fftscore );
#endif
	if( localalloclen < nlen )
	{
		if( localalloclen )
		{
			if( !kobetsubunkatsu )
			{
				FreeFukusosuuMtx ( seqVector1 );
				FreeFukusosuuMtx ( seqVector2 );
				FreeFukusosuuVec( naisekiNoWa );
				FreeFukusosuuMtx( naiseki );
				FreeDoubleVec( soukan );
			}
		}

		if( !kobetsubunkatsu )
		{
			naisekiNoWa = AllocateFukusosuuVec( nlen );
			naiseki = AllocateFukusosuuMtx( n20or4or2, nlen );
			seqVector1 = AllocateFukusosuuMtx( n20or4or2+1, nlen+1 );
			seqVector2 = AllocateFukusosuuMtx( n20or4or2+1, nlen+1 );
			soukan = AllocateDoubleVec( nlen+1 );
		}
		localalloclen = nlen;
	}

	for( j=0; j<clus1; j++ ) strcpy( tmpseq1[j], seq1[j] );
	for( j=0; j<clus2; j++ ) strcpy( tmpseq2[j], seq2[j] );

#if DEBUG
fftfp = fopen( "input_of_Falign", "w" );
fprintf( fftfp, "nlen = %d\n", nlen );
fprintf( fftfp, "seq1: ( %d sequences ) \n", clus1 );
for( i=0; i<clus1; i++ )
	fprintf( fftfp, "%s\n", seq1[i] );
fprintf( fftfp, "seq2: ( %d sequences ) \n", clus2 );
for( i=0; i<clus2; i++ )
	fprintf( fftfp, "%s\n", seq2[i] );
fclose( fftfp );
system( "less input_of_Falign < /dev/tty > /dev/tty" );
#endif
	if( !kobetsubunkatsu )
	{
		if( fftkeika ) fprintf( stderr,  " FFT ... " );

		for( j=0; j<n20or4or2; j++ ) vec_init( seqVector1[j], nlen );
		if( fftscore && scoremtx != -1 )
		{
			for( i=0; i<clus1; i++ )
			{
				seq_vec_5( seqVector1[0], polarity, volume, eff1[i], tmpseq1[i] );
			}
		}
		else
		{
			for( i=0; i<clus1; i++ )
				seq_vec_3( seqVector1, eff1[i], tmpseq1[i] );
		}
#if DEBUG
fftfp = fopen( "seqVec", "w" );
fprintf( fftfp, "before transform\n" );
for( k=0; k<n20or4or2; k++ ) 
{
   fprintf( fftfp, "nlen=%d\n", nlen );
   fprintf( fftfp, "%c\n", amino[k] );
   for( l=0; l<nlen; l++ )
   fprintf( fftfp, "%f %f\n", seqVector1[k][l].R, seqVector1[k][l].I );
}
fprintf(fftfp, "sequences weight: \n");
for(k = 0; k < clus1; ++ k) fprintf(fftfp, "%f ", eff1[k]);
fprintf(fftfp, "\n");
fclose( fftfp );
//system( "less seqVec < /dev/tty > /dev/tty" );
#endif

		for( j=0; j<n20or4or2; j++ ) vec_init( seqVector2[j], nlen );
		if( fftscore && scoremtx != -1 )
		{
			for( i=0; i<clus2; i++ )
			{
				seq_vec_5( seqVector2[0], polarity, volume, eff2[i], tmpseq2[i] );
			}
		}
		else
		{
			for( i=0; i<clus2; i++ )
				seq_vec_3( seqVector2, eff2[i], tmpseq2[i] );
		}

#if DEBUG
fftfp = fopen( "seqVec2", "w" );
fprintf( fftfp, "before fft\n" );
for( k=0; k<n20or4or2; k++ ) 
{
   fprintf( fftfp, "%c\n", amino[k] );
   for( l=0; l<nlen; l++ )
   fprintf( fftfp, "%f %f\n", seqVector2[k][l].R, seqVector2[k][l].I );
}
fprintf(fftfp, "sequences weight: \n");
for(k = 0; k < clus2; ++ k) fprintf(fftfp, "%f ", eff2[k]);
fprintf(fftfp, "\n");
fclose( fftfp );
//system( "less seqVec2 < /dev/tty > /dev/tty" );
#endif

		for( j=0; j<n20or4or2; j++ )
		{
			fft( nlen, seqVector2[j], 0 );
			fft( nlen, seqVector1[j], 0 );
		}
#if DEBUG
fftfp = fopen( "seqVec3", "w" );
fprintf( fftfp, "after fft\n" );
for( k=0; k<n20or4or2; k++ ) 
{
   fprintf( fftfp, "%c\n", amino[k] );
   for( l=0; l<nlen; l++ )
	   fprintf( fftfp, "%f %f\n", seqVector2[k][l].R, seqVector2[k][l].I );
}
fclose( fftfp );
//system( "less seqVec2 < /dev/tty > /dev/tty" );
#endif

		for( k=0; k<n20or4or2; k++ ) 
		{
			for( l=0; l<nlen; l++ ) 
				calcNaiseki( naiseki[k]+l, seqVector1[k]+l, seqVector2[k]+l );
		}
		for( l=0; l<nlen; l++ ) 
		{
			naisekiNoWa[l].R = 0.0;
			naisekiNoWa[l].I = 0.0;
			for( k=0; k<n20or4or2; k++ ) 
			{
				naisekiNoWa[l].R += naiseki[k][l].R;
				naisekiNoWa[l].I += naiseki[k][l].I;
			}
		}
	
#if DEBUG
	fftfp = fopen( "naisekiNoWa", "w" );
	fprintf( fftfp, "#Before fft\n" );
	for( l=0; l<nlen; l++ )
		fprintf( fftfp, "%d  %f %f\n", l, naisekiNoWa[l].R, naisekiNoWa[l].I ); 
	fclose( fftfp );
	//system( "less naisekiNoWa < /dev/tty > /dev/tty " );
#endif

		fft( -nlen, naisekiNoWa, 0 );
	
		for( m=0; m<=nlen2; m++ ) 
			soukan[m] = naisekiNoWa[nlen2-m].R;
		for( m=nlen2+1; m<nlen; m++ ) 
			soukan[m] = naisekiNoWa[nlen+nlen2-m].R;

#if DEBUG
	fftfp = fopen( "naisekiNoWa", "a" );
	fprintf( fftfp, "#After fft\n" );
	for( l=0; l<nlen; l++ )
		fprintf( fftfp, "%d  %f\n", l, naisekiNoWa[l].R ); 
	fclose( fftfp );
	fftfp = fopen( "list.plot", "w"  );
	fprintf( fftfp, "plot 'naisekiNoWa'\npause -1" );
	fclose( fftfp );
	//system( "/usr/bin/gnuplot list.plot &" );
#endif
#if DEBUG
	fftfp = fopen("soukan", "w");
	fprintf( fftfp, "soukan\n" );
	for( l=0; l<nlen; l++ )
		fprintf( fftfp, "%d  %f\n", l-nlen2, soukan[l] ); 
#if DEBUG
	fftfp = fopen( "list.plot", "w"  );
	fprintf( fftfp, "plot 'frt'\n pause +1" );
	fclose( fftfp );
	//system( "/usr/bin/gnuplot list.plot" );
#endif
#endif
		getKouho( kouho, NKOUHO, soukan, nlen );
	}
	count = 0;
	if( kobetsubunkatsu )
	{
		maxk = 1;
		kouho[0] = 0;
	}
	else
	{
		maxk = NKOUHO;
	}
	for( k=0; k<maxk; k++ ) 
	{
		lag = kouho[k];
		if( lag <= -len1 || len2 <= lag ) continue;
		zurasu2( lag, clus1, clus2, seq1, seq2, tmpptr1, tmpptr2 );
		tmpint = alignableReagion( clus1, clus2, tmpptr1, tmpptr2, eff1, eff2, segment+count );
		if( count+tmpint > MAXSEG -3 ) ErrorExit( "TOO MANY SEGMENTS.\n" );


		if( tmpint == 0 ) break; // 060430 iinoka ?
		while( tmpint-- > 0 )
		{
			if( lag > 0 )
			{
				segment1[count].start  = segment[count].start ;
				segment1[count].end    = segment[count].end   ;
				segment1[count].center = segment[count].center;
				segment1[count].score  = segment[count].score;

				segment2[count].start  = segment[count].start  + lag;
				segment2[count].end    = segment[count].end    + lag;
				segment2[count].center = segment[count].center + lag;
				segment2[count].score  = segment[count].score       ;
			}
			else
			{
				segment1[count].start  = segment[count].start  - lag;
				segment1[count].end    = segment[count].end    - lag;
				segment1[count].center = segment[count].center - lag;
				segment1[count].score  = segment[count].score       ;

				segment2[count].start  = segment[count].start ;
				segment2[count].end    = segment[count].end   ;
				segment2[count].center = segment[count].center;
				segment2[count].score  = segment[count].score ;
			}
#if DEBUG
			fprintf( stderr, "Segment: \n");
			fprintf( stderr, "in 1 %d\n", segment1[count].center );
			fprintf( stderr, "in 2 %d\n", segment2[count].center );
#endif
			segment1[count].pair = &segment2[count];
			segment2[count].pair = &segment1[count];
			count++;
		}
	}
#if DEBUG
	if( !kobetsubunkatsu && fftkeika )
		fprintf( stderr, "%d anchors found\r", count );
#endif
	if( !count && fftNoAnchStop )
		ErrorExit( "Cannot detect anchor!" );
#if DEBUG
	fftfp = fopen( "homo", "w" );
	fprintf( fftfp, "RESULT before sort:\n" );
	for( l=0; l<count+1; l++ )
	{
		fprintf( fftfp, "cut[%d]=%d, ", l, segment1[l].center );
		fprintf( fftfp, "%d score = %f\n", segment2[l].center, segment1[l].score );
	}
	fclose(fftfp);
#endif

#if KEIKA
	fprintf( stderr, "done. (%d anchors)\n", count );
	fprintf( stderr, "Aligning anchors ... " );
#endif
	for( i=0; i<count; i++ )
	{
		sortedseg1[i] = &segment1[i];
		sortedseg2[i] = &segment2[i];
	}
	mymergesort( 0, count-1, sortedseg1 ); 
	mymergesort( 0, count-1, sortedseg2 ); 
	for( i=0; i<count; i++ ) sortedseg1[i]->number = i;
	for( i=0; i<count; i++ ) sortedseg2[i]->number = i;


	if( kobetsubunkatsu )
	{
		for( i=0; i<count; i++ )
		{
			cut1[i+1] = sortedseg1[i]->center;
			cut2[i+1] = sortedseg2[i]->center;
		}
		cut1[0] = 0;
		cut2[0] = 0;
		cut1[count+1] = len1;
		cut2[count+1] = len2;
		count += 2;
	}
	else
	{
		if( crossscoresize < count+2 )
		{
			crossscoresize = count+2;
#if 1
			if( fftkeika ) fprintf( stderr, "######allocating crossscore, size = %d\n", crossscoresize );
#endif
			if( crossscore ) FreeDoubleMtx( crossscore );
			crossscore = AllocateDoubleMtx( crossscoresize, crossscoresize );
		}
		for( i=0; i<count+2; i++ ) for( j=0; j<count+2; j++ )
			crossscore[i][j] = 0.0;
		for( i=0; i<count; i++ )
		{
			crossscore[segment1[i].number+1][segment1[i].pair->number+1] = segment1[i].score;
			cut1[i+1] = sortedseg1[i]->center;
			cut2[i+1] = sortedseg2[i]->center;
		}

#if DEBUG
		fftfp = fopen( "homo", "a" );
		fprintf( fftfp, "AFTER SORT\n" );
		for( i=0; i<count+1; i++ ) fprintf( fftfp, "%d, %d\n", cut1[i], cut2[i] );
		fclose(fftfp);
		//fprintf( stderr, "crossscore = \n" );
		//for( i=0; i<count+1; i++ )
		//{
		//	for( j=0; j<count+1; j++ )
		//		fprintf( stderr, "%.0f ", crossscore[i][j] );
		//	fprintf( stderr, "\n" );
		//}
#endif

		crossscore[0][0] = 10000000.0;
		cut1[0] = 0; 
		cut2[0] = 0;
		crossscore[count+1][count+1] = 10000000.0;
		cut1[count+1] = len1;
		cut2[count+1] = len2;
		count += 2;
		count0 = count;
	
		blockAlign2( cut1, cut2, sortedseg1, sortedseg2, crossscore, &count );

//		if( count-count0 )
//			fprintf( stderr, "%d unused anchors\n", count0-count );

		if( !kobetsubunkatsu && fftkeika )
			fprintf( stderr, "%d anchors found\n", count );
		if( fftkeika )
		{
			if( count0 > count )
			{
				fprintf( stderr, "REPEAT!? \n" ); 
				if( fftRepeatStop ) exit( 1 );
			}
#if KEIKA
			else fprintf( stderr, "done\n" );
#endif
		}
	}

#if DEBUG
	fftfp = fopen( "homo", "a" );
	fprintf( fftfp, "RESULT after sort:\n" );
	for( l=0; l<count; l++ )
	{
		fprintf( fftfp, "cut[%d]=%d, ", l, segment1[l].center );
		fprintf( fftfp, "%d\n", segment2[l].center );
	}
	fclose( fftfp );
#endif

#if DEBUG
	fprintf( stderr, "RESULT after blockalign:\n" );
	for( l=0; l<count+1; l++ )
	{
		fprintf( stderr, "cut : %d %d\n", cut1[l], cut2[l] );
	}
#endif

#if DEBUG
	fprintf( trap_g, "Divided to %d segments\n", count-1 );
	fprintf( trap_g, "%d %d %d forg\n",  clus1, clus2 , count-1 );
#endif

	totallen = 0;
	for( j=0; j<clus1; j++ ) result1[j][0] = 0;
	for( j=0; j<clus2; j++ ) result2[j][0] = 0;
	totalscore = 0.0;
	*fftlog = -1;
#if DEBUG
	reporterr( "\nin Falign(), *fftlog = %d\n", *fftlog );
#endif
	for( i=0; i<count-1; i++ )
	{
		*fftlog += 1;
		if( i == 0 ) headgp = outgap; else headgp = 1;
		if( i == count-2 ) tailgp = outgap; else tailgp = 1;

		if( cut1[i] )
			getkyokaigap( sgap1, seq1, cut1[i]-1, clus1 );
		else
			for( j=0; j<clus1; j++ ) sgap1[j] = 'o';
		if( cut2[i] )
			getkyokaigap( sgap2, seq2, cut2[i]-1, clus2 );
		else
			for( j=0; j<clus2; j++ ) sgap2[j] = 'o';

		if( cut1[i+1] != len1 )
			getkyokaigap( egap1, seq1, cut1[i+1], clus1 );
		else    
			for( j=0; j<clus1; j++ ) egap1[j] = 'o';
		if( cut2[i+1] != len2 )
			getkyokaigap( egap2, seq2, cut2[i+1], clus2 );
		else    
			for( j=0; j<clus2; j++ ) egap2[j] = 'o';
#if DEBUG
		{
			fprintf( stderr, "kyokkaigap1(%d)=", cut1[i]-1 );
			for( j=0; j<clus1; j++ )
				fprintf( stderr, "%c", sgap1[j] );
			fprintf( stderr, "=kyokkaigap1-start\n" );
		}
		{
			fprintf( stderr, "kyokkaigap2(%d)=", cut2[i]-1 );
			for( j=0; j<clus2; j++ )
				fprintf( stderr, "%c", sgap2[j] );
			fprintf( stderr, "=kyokkaigap2-start\n" );
		}
		{
			fprintf( stderr, "kyokkaigap1(%d)=", cut1[i]-1 );
			for( j=0; j<clus1; j++ )
				fprintf( stderr, "%c", egap1[j] );
			fprintf( stderr, "=kyokkaigap1-end\n" );
		}
		{
			fprintf( stderr, "kyokkaigap2(%d)=", cut2[i]-1 );
			for( j=0; j<clus2; j++ )
				fprintf( stderr, "%c", egap2[j] );
			fprintf( stderr, "=kyokkaigap2-end\n" );
		}
#endif

#if DEBUG
		fprintf( stderr, "DP %03d / %03d %4d to ", i+1, count-1, totallen );
#endif

		// reporterr( "cut1[] = %d\n", cut1[i] );
		// reporterr( "cut2[] = %d\n", cut2[i] );

		for( j=0; j<clus1; j++ )
		{
			strncpy( tmpres1[j], seq1[j]+cut1[i], cut1[i+1]-cut1[i] );
			tmpres1[j][cut1[i+1]-cut1[i]] = 0;
		}
		if( kobetsubunkatsu && fftkeika ) commongappick( clus1, tmpres1 ); //dvtditr に呼ばれたとき fftkeika=1

		for( j=0; j<clus2; j++ )
		{
			strncpy( tmpres2[j], seq2[j]+cut2[i], cut2[i+1]-cut2[i] );
			tmpres2[j][cut2[i+1]-cut2[i]] = 0;
		}
		if( kobetsubunkatsu && fftkeika ) commongappick( clus2, tmpres2 ); //dvtditr に呼ばれたとき fftkeika=1

		if( constraint )
		{
			fprintf( stderr, "Not supported\n" );
			exit( 1 );
		}
#if 0
		fprintf( stderr, "i=%d, before alignment", i );
		fprintf( stderr, "%4d\n", totallen );
		fprintf( stderr, "\n\n" );
		for( j=0; j<clus1; j++ ) 
		{
			fprintf( stderr, "%s\n", tmpres1[j] );
		}
		fprintf( stderr, "-------\n" );
		for( j=0; j<clus2; j++ ) 
		{
			fprintf( stderr, "%s\n", tmpres2[j] );
		}
#endif

#if DEBUG
		fprintf( stdout, "writing input\n" );
		for( j=0; j<clus1; j++ )
		{
			fprintf( stdout, ">%d of GROUP1\n", j );
			fprintf( stdout, "%s\n", tmpres1[j] );
		}
		for( j=0; j<clus2; j++ )
		{
			fprintf( stdout, ">%d of GROUP2\n", j );
			fprintf( stdout, "%s\n", tmpres2[j] );
		}
		fflush( stdout );
#endif
		switch( alg )
		{
			case( 'a' ):
				totalscore += Aalign( tmpres1, tmpres2, eff1, eff2, clus1, clus2, alloclen );
				break;
			case( 'A' ): // Normal mode
				if( clus1 == 1 && clus2 == 1 ) // two sequences alignment
					totalscore += G__align11( n_dynamicmtx, tmpres1, tmpres2, alloclen, headgp, tailgp );
				else // Multiple sequence alignment
					totalscore += A__align( n_dynamicmtx, penalty, penalty_ex, tmpres1, tmpres2, eff1, eff2, clus1, clus2, alloclen, sgap1, sgap2, egap1, egap2, headgp, tailgp);
				break;
			default:
				fprintf( stderr, "alg = %c\n", alg );
				ErrorExit( "NOT SUPPORTED FEATURE IN Falign.c" );
				break;
		}

#ifdef enablemultithread
		if( chudanres && *chudanres )
		{
//			fprintf( stderr, "\n\n## CHUUDAN!!! at Falign_localhom\n" );
			return( -1.0 );
		}
#endif

		nlen = strlen( tmpres1[0] );
		if( totallen + nlen > alloclen )
		{
			fprintf( stderr, "totallen=%d +  nlen=%d > alloclen = %d\n", totallen, nlen, alloclen );
			ErrorExit( "LENGTH OVER in Falign\n " );
		}
		for( j=0; j<clus1; j++ ) strcat( result1[j], tmpres1[j] );
		for( j=0; j<clus2; j++ ) strcat( result2[j], tmpres2[j] );
		totallen += nlen;
#if 0
		fprintf( stderr, "$#####$$$$ i=%d", i );
		fprintf( stderr, "%4d\n", totallen );
		fprintf( stderr, "\n\n" );
		for( j=0; j<clus1; j++ ) 
		{
			fprintf( stderr, "%s\n", tmpres1[j] );
		}
		fprintf( stderr, "-------\n" );
		for( j=0; j<clus2; j++ ) 
		{
			fprintf( stderr, "%s\n", tmpres2[j] );
		}
#endif
	}
	//reporterr( "\nafter Falign(), *fftlog = %d\n", *fftlog );

#if KEIKA
	fprintf( stderr, "DP ... done   \n" );
#endif

	for( j=0; j<clus1; j++ ) strcpy( seq1[j], result1[j] );
	for( j=0; j<clus2; j++ ) strcpy( seq2[j], result2[j] );
#if DEBUG
	for( j=0; j<clus1; j++ ) 
	{
		fprintf( stderr, "in Falign, %s\n", result1[j] );
	}
	fprintf( stderr, "- - - - - - - - - - -\n" );
	for( j=0; j<clus2; j++ ) 
	{
		fprintf( stderr, "in Falign, %s\n", result2[j] );
	}
#endif


	FreeCharMtx( result1 );
	FreeCharMtx( result2 );
	FreeCharMtx( tmpres1 );
	FreeCharMtx( tmpres2 );
	free( sgap1 );
	free( egap1 );
	free( sgap2 );
	free( egap2 );
	FreeCharMtx( tmpseq1 );
	FreeCharMtx( tmpseq2 );
	free( tmpptr1 );
	free( tmpptr2 );
	return( totalscore );
}
