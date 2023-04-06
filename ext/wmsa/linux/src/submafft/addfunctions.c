#include "mltaln.h"

static void strncpy0( char *s1, char *s2, int n )
{
	while( n-- ) *s1++ = *s2++;
	*s1 = 0;
}

#if 0
static void strncpy0x( char *s1, char *s2, int n ) { while( n-- ) *s1++ = *s2++; *s1 = 0; }
static void strncpy0b0( char *s1, char *s2, int n ) { while( n-- ) *s1++ = *s2++; *s1 = 0; }
static void strncpy0b1( char *s1, char *s2, int n ) { while( n-- ) *s1++ = *s2++; *s1 = 0; }
static void strncpy0b2( char *s1, char *s2, int n ) { while( n-- ) *s1++ = *s2++; *s1 = 0; }
static void strncpy0n0( char *s1, char *s2, int n ) { while( n-- ) *s1++ = *s2++; *s1 = 0; }
static void strncpy0n1( char *s1, char *s2, int n ) { while( n-- ) *s1++ = *s2++; *s1 = 0; }
static void strncpy0n2( char *s1, char *s2, int n ) { while( n-- ) *s1++ = *s2++; *s1 = 0; }
static void strncpy0a0( char *s1, char *s2, int n ) { while( n-- ) *s1++ = *s2++; *s1 = 0; }
static void strncpy0a1( char *s1, char *s2, int n ) { while( n-- ) *s1++ = *s2++; *s1 = 0; }
static void strncpy0a2( char *s1, char *s2, int n ) { while( n-- ) *s1++ = *s2++; *s1 = 0; }
static void strncpy0o0( char *s1, char *s2, int n ) { while( n-- ) *s1++ = *s2++; *s1 = 0; }
static void strncpy0o1( char *s1, char *s2, int n ) { while( n-- ) *s1++ = *s2++; *s1 = 0; }
static void strncpy0o2( char *s1, char *s2, int n ) { while( n-- ) *s1++ = *s2++; *s1 = 0; }
#endif

static void eqpick( char *aseq, char *seq )
{
	for( ; *seq != 0; seq++ )
	{
		if( *seq != '=' )
			*aseq++ = *seq;
	}
	*aseq = 0;

}


void eq2dashmatomete( char **s, int n )
{
	int i, j;
	char sj;

	for( j=0; (sj=s[0][j]); j++ )
	{
		if( sj == '=' )
		{
			for( i=0; i<n; i++ )
			{
				s[i][j] = '-';
			}
		}
	}
}

void eq2dashmatometehayaku( char **s, int n )
{
	int i, j, c;
	int *tobechanged;
	int len = strlen( s[0] );

	tobechanged = calloc( len+1, sizeof( int ) ); // len+1, 2017/Nov/15
	c = 0;
	for( j=0; j<len; j++ )
	{
		if( s[0][j] == '=' ) tobechanged[c++] = j;
	}
	tobechanged[c] = -1;

	for( i=0; i<n; i++ )
	{
		for( c=0; (j=tobechanged[c])!=-1; c++ )
			s[i][j] = '-';
	}
	free( tobechanged );
}

void eq2dash( char *s )
{
	while( *s )
	{
		if( *s == '=' ) 
		{
			*s = '-';
		}
		s++;
	}
}

static void plus2gapchar( char *s, char gapchar )
{
	while( *s )
	{
		if( *s == '+' ) 
		{
			*s = gapchar;
		}
		s++;
	}
}

void findnewgaps( int n, int rep, char **seq, int *gaplen )
{
	int i, pos, len, len1;

	len = strlen( seq[0] );	
//	for( i=0; i<len; i++ ) gaplen[i] = 0; // calloc de shokika sareteirukara hontou ha iranai
	len1 = len + 1;
	for( i=0; i<len1; i++ ) gaplen[i] = 0; // realloc de shokika sareteirukara iru!
	
	pos = 0;
	for( i=0; i<len; i++ )
	{
		if( seq[rep][i] == '=' ) 
		{
			if( disp ) fprintf( stderr, "Newgap! pos = %d\n", pos );
			gaplen[pos]++;
		}
		else
			pos++;
	}

#if 0
	if( disp )
	{
		fprintf( stderr, "\ngaplen[] in findnewgaps() = \n" );
		for(i=0; i<pos; i++ ) fprintf( stderr, "%d ", gaplen[i] );
		fprintf( stderr, "\n" );
		reporterr( "pos=%d\n", pos );
	}
#endif
}

void findcommongaps( int n, char **seq, int *gapmap )
{
	int i, j, pos, len, len1;
	len = strlen( seq[0] );	
	len1 = len+1;

//	fprintf( stderr, "seq[0] = %s\n", seq[0] );
	for( i=0; i<len1; i++ ) gapmap[i] = 0;
	
	pos = 0;
	for( i=0; i<len; i++ )
	{
		for( j=0; j<n; j++ )
			if( seq[j][i] != '-' ) break;

		if( j == n ) gapmap[pos]++;
		else
			pos++;
	}
#if 0
	for( i=0; i<pos; i++ )
	{
		fprintf( stderr, "vec[%d] = %d\n", i, gapmap[i] );
	}
#endif
}

void adjustgapmap( int newlen, int *gapmap, char *seq )
{
	int j;
	int pos;
	int newlen1 = newlen+1;
	int *tmpmap;


	tmpmap = AllocateIntVec( newlen+2 );
	j = 0;
	pos = 0;
	while( *seq )
	{
//		fprintf( stderr, "j=%d *seq = %c\n", j, *seq );
		if( *seq++ == '=' )
			tmpmap[j++] = 0;
		else
		{
			tmpmap[j++] = gapmap[pos++];
		}
	}
	tmpmap[j++] = gapmap[pos];

	for(j=0; j<newlen1; j++)
		gapmap[j] = tmpmap[j];

	free( tmpmap );

#if 0
	reporterr( "gapmap in adjustgapmap() =\n" );
	for(j=0; j<newlen1; j++) reporterr( "%d ", gapmap[j] );
	reporterr( "length = %d\n", newlen );
#endif
}


static int countnogaplen( int *gaplen, int *term )
{
	int v = 0;
	while( gaplen < term )
	{
		if( *gaplen++ == 0 ) v++;
		else break;
	}
	return( v );
}

static int countgapmap( int *gapmap, int *term )
{
	int v = 0;
	while( gapmap < term )
	{
//		reporterr( "*gapmap = %d\n", *gapmap );
		if( *gapmap++ == 0 ) v++;
		else break;
	}
	return( v );
}



static void reflectsmoothing( char *ref, int *mem, char **seq, int len )
{
	char *tmpseq;
	int i, j, k, p;

//	reporterr( "#### reflectsmoothing!!!!!\n" );

//	if( mem[1] != -1 ) reporterr( "original = %s\n", seq[mem[1]] );

	tmpseq = calloc( len+1, sizeof( char ) );


	for( j=1; (i=mem[j])!=-1; j++ )
	{
		eqpick( tmpseq, seq[i] );
		for( k=0, p=0; p<len; k++ )
		{
			while( ref[p] == '=' ) seq[i][p++] = '=';
			seq[i][p++] = tmpseq[k];
		}
	}
	free( tmpseq );

//	if( mem[1] != -1 ) reporterr( "output   = %s\n", seq[mem[1]] );

//	reporterr( "#### done!!!!!\n" );
}

static int smoothing1rightmulti( int len, char *ref ) // osoi!
{
	int i, j, k;
	int shiftfrom = -1;
	int shiftto = -1;
	int *hit;
	int val = 0, nhit = 0;

	hit = NULL;

//	reporterr( "ref (1rightmulti) = %s\n", ref );

	for( i=1, nhit=0; i<len-1; i++ ) // break nashi no baai, hidarihaji ha saigo
//	for( i=len-2; i>0; i-- ) // break ari no baai, migihajiha saigo
	{
		if( ref[i-1] == '+' && ( ref[i] != '+' && ref[i] != '=' ) && ref[i+1] == '=' )
		{
//			reporterr( "hit! i=%d, len=%d\n", i, len );
			hit = realloc( hit, (nhit+1) * sizeof( int ) );
			hit[nhit] = i;
			nhit += 1;
//			break;
		}
	}
	if( nhit == 0 ) return( 0 );


	for( k=0; k<nhit; k++ )
	{
		for( j=hit[k]+1; j<=len; j++ )
		{
			if( ref[j] != '=' )
			{
				shiftto = j-1;
				break;
			}
		}
		if( j == len && ref[len-1] == '=' )
		{
			reporterr( "hit[i].end = %d, j = len-1, skip!\n" );
			continue;
		}

		if( shiftto < len-1 && ref[shiftto+1] == '+' ) continue; // muda dakara

		val += 1;
		shiftfrom = hit[k];
		if( ref[shiftto] != '=' ) // atode sakujo 
		{
			reporterr( "Error in smoothing1left!\n" );
			exit( 1 );
		}
		ref[shiftto] = ref[shiftfrom];
		ref[shiftfrom] = '=';
	}
	free( hit );

//	reporterr( "ref (1rightmulti) = %s\n", ref );
	reporterr( " %d out of %d have been smoothed (right).\n", val, nhit );

//	if( nhit > 1 ) exit( 1 );
	return( val );
}

static int smoothing1leftmulti( int len, char *ref ) // osoi!
{
	int i, j, k;
	int shiftfrom = -1;
	int shiftto = -1;
	int *hit;
	int val = 0, nhit = 0;

	hit = NULL;

//	reporterr( "ref (1leftmulti) = %s\n", ref );

	for( i=1, nhit=0; i<len-1; i++ ) // break nashi no baai, hidarihaji ha saigo
//	for( i=len-2; i>0; i-- ) // break ari no baai, migihajiha saigo
	{
		if( ref[i-1] == '=' && ( ref[i] != '+' && ref[i] != '=' ) && ref[i+1] == '+' )
		{
//			reporterr( "hit! i=%d, len=%d\n", i, len );
			hit = realloc( hit, (nhit+1) * sizeof( int ) );
			hit[nhit] = i;
			nhit += 1;
//			break;
		}
	}
	if( nhit == 0 ) return( 0 );

	for( k=0; k<nhit; k++ )
	{
		for( j=hit[k]-1; j>-1; j-- )
		{
			if( ref[j] != '=' )
			{
				shiftto = j+1;
				break;
			}
		}
		if( j == -1 && ref[0] == '=' )
		{
			reporterr( "hit[i].end = %d, j = -1, skip!\n" );
			continue;
		}

		if( shiftto > 0 && ref[shiftto-1] == '+' ) continue; // muda dakara

		val += 1;
		shiftfrom = hit[k];
		if( ref[shiftto] != '=' ) // atode sakujo 
		{
			reporterr( "Error in smoothing1left!\n" );
			exit( 1 );
		}
		ref[shiftto] = ref[shiftfrom];
		ref[shiftfrom] = '=';
	
	}
	free( hit );

//	reporterr( "ref (1leftmulti) = %s\n", ref );
	reporterr( " %d out of %d have been smoothed (left).\n", val, nhit );


//	if( nhit > 1 ) exit( 1 );
	return( val );
}

void restorecommongapssmoothly( int njob, int n0, char **seq, int *ex1, int *ex2, int *gapmap, int alloclen, char gapchar )
{
	int *mem;
	char *tmpseq;
	char *cptr;
	int *iptr;
	int *tmpgapmap;
	int i, j, k, len, rep1, rep2, len1, klim, leninserted;
	int totalres;
	
	if( n0 == 0 ) return;


	mem = calloc( njob+1, sizeof( int ) ); // +1 ha iranai.
	intcpy( mem, ex1 );
	intcat( mem, ex2 );
//	tmpseq = calloc( alloclen+2, sizeof( char ) );
//	tmpgapmap = calloc( alloclen+2, sizeof( int ) );

#if 0 // iranai
	for( i=0; (k=mem[i])!=-1; i++ ) // iranai
		reporterr( "mem[%d] = %d\n", i, k ); // iranai
	if( i == njob ) // iranai
	{
		fprintf( stderr, "Error in restorecommongaps()\n" );
		free( mem );
		exit( 1 );
	}
#endif
	rep1 = ex1[0];
	rep2 = ex2[0];
	len = strlen( seq[rep1] );
	len1 = len+1;

	tmpseq = calloc( alloclen, sizeof( char ) );
	tmpgapmap = calloc( alloclen, sizeof( int ) );

#if 0
	reporterr( "\n" );
	reporterr( "seq[rep1] = %s\n", seq[rep1] );
	reporterr( "seq[rep2] = %s\n", seq[rep2] );
#endif

	for( k=0; (i=mem[k])!=-1; k++ )
	{
		cptr = tmpseq;
		for( j=0; j<len1; j++ )
		{
			klim = gapmap[j];
//			for( k=0; k<gapmap[j]; k++ )
			while( klim-- )
				*(cptr++) = '+'; // ???
			*(cptr++) = seq[i][j];
		}
		*cptr = 0;
		strcpy( seq[i], tmpseq );
	}
#if 0
	reporterr( "->\n" );
	reporterr( "seq[rep1] = \n%s\n", seq[rep1] );
	reporterr( "seq[rep2] = \n%s\n", seq[rep2] );
#endif

	leninserted = strlen( seq[rep1] );
#if 0
	reporterr( "gapmap =\n" );
	for(j=0; j<len1; j++) 
	{
		reporterr( "%d", gapmap[j] );
		for( i=gapmap[j]; i>0; i-- ) reporterr( "-" );
	}
	reporterr( "\n" );
#endif

#if 0
	resprev = 10000; // tekitou
	while( 1 )
	{
		res = 0;
//		reporterr( "\nsmoothing1right..\n" );
		res  = (0<smoothing1right( leninserted, seq[rep1], gapmap, seq, ex1 ));
//		reporterr( "done. res = %d\n", res );
//		reporterr( "smoothing1right..\n" );
		res += (0<smoothing1right( leninserted, seq[rep2], gapmap, seq, ex2 ));
//		reporterr( "done. res = %d\n", res );

//		reporterr( "smoothing1left..\n" );
		res += (0<smoothing1left( leninserted, seq[rep1], gapmap, seq, ex1 ));
//		reporterr( "done. res = %d\n", res );
//		reporterr( "smoothing1left..\n" );
		res += (0<smoothing1left( leninserted, seq[rep2], gapmap, seq, ex2 ));
//		reporterr( "done. res = %d\n", res );

		reporterr( " Smoothing .. %d \n", res );
		if( res >= resprev ) break;
//		if( res == 0 ) break;
		resprev = res;
	}
#else
	totalres = 0;
	totalres += smoothing1rightmulti( leninserted, seq[rep1] );
	totalres += smoothing1leftmulti( leninserted, seq[rep1] );
	if( totalres ) reflectsmoothing( seq[rep1], ex1, seq, leninserted );

	totalres = 0;
	totalres += smoothing1rightmulti( leninserted, seq[rep2] );
	totalres += smoothing1leftmulti( leninserted, seq[rep2] );
	if( totalres ) reflectsmoothing( seq[rep2], ex2, seq, leninserted );
#endif

	for( k=0; (i=mem[k])!=-1; k++ ) plus2gapchar( seq[i], gapchar );

#if 0
	reporterr( "->\n" );
	reporterr( "seq[rep1] = \n%s\n", seq[rep1] );
	reporterr( "seq[rep2] = \n%s\n", seq[rep2] );
	reporterr( "gapmap =\n" );
	for(j=0; j<len1; j++) 
	{
		reporterr( "%d", gapmap[j] );
		for( i=gapmap[j]; i>0; i-- ) reporterr( "-" );
	}
	reporterr( "\n" );
#endif

	iptr = tmpgapmap;
	for( j=0; j<len1; j++ )
	{
		*(iptr++) = gapmap[j];
		for( k=0; k<gapmap[j]; k++ )
			*(iptr++) = 0;
	}
	*iptr = -1;

	intcpy( gapmap, tmpgapmap );
//	iptr = tmpgapmap;
//	while( *iptr != -1 ) *gapmap++ = *iptr++;

	free( mem );
	free( tmpseq );
	free( tmpgapmap );
}

void restorecommongaps( int njob, int n0, char **seq, int *ex1, int *ex2, int *gapmap, int alloclen, char gapchar )
{
	int *mem;
	char *tmpseq;
	char *cptr;
	int *iptr;
	int *tmpgapmap;
	int i, j, k, len, rep, len1, klim;
	

	if( n0 == 0 ) return;


	mem = calloc( njob+1, sizeof( int ) ); // +1 ha iranai.
	intcpy( mem, ex1 );
	intcat( mem, ex2 );
//	tmpseq = calloc( alloclen+2, sizeof( char ) );
//	tmpgapmap = calloc( alloclen+2, sizeof( int ) );

#if 0 // iranai
	for( i=0; (k=mem[i])!=-1; i++ ) // iranai
		reporterr( "mem[%d] = %d\n", i, k ); // iranai
	if( i == njob ) // iranai
	{
		fprintf( stderr, "Error in restorecommongaps()\n" );
		free( mem );
		exit( 1 );
	}
#endif
	rep = mem[0];
	len = strlen( seq[rep] );
	len1 = len+1;

	tmpseq = calloc( alloclen, sizeof( char ) );
	tmpgapmap = calloc( alloclen, sizeof( int ) );



	for( k=0; (i=mem[k])!=-1; k++ )
	{
		cptr = tmpseq;
		for( j=0; j<len1; j++ )
		{
			klim = gapmap[j];
//			for( k=0; k<gapmap[j]; k++ )
			while( klim-- )
				*(cptr++) = gapchar; // ???
			*(cptr++) = seq[i][j];
		}
		*cptr = 0;
		strcpy( seq[i], tmpseq );
	}

	iptr = tmpgapmap;
	for( j=0; j<len1; j++ )
	{
		*(iptr++) = gapmap[j];
		for( k=0; k<gapmap[j]; k++ )
			*(iptr++) = 0;
	}
	*iptr = -1;

	iptr = tmpgapmap;
	while( *iptr != -1 ) *gapmap++ = *iptr++;

	free( mem );
	free( tmpseq );
	free( tmpgapmap );
}

int deletenewinsertions_whole_eq( int on, int an, char **oseq, char **aseq, int **deletelist )
{
	int i, j, p, q, allgap, ndel;
	int len = strlen( oseq[0] );
	char *eqseq, tmpc;

//	reporterr( "In deletenewinsertions_whole_eq\n" );
//	for( j=0; j<on; j++ ) reporterr( "\no=%s\n", oseq[j] );
//	for( j=0; j<an; j++ ) reporterr( "a=%s\n", aseq[j] );

	eqseq = calloc( len+1, sizeof( char ) );
	for( i=0; i<len; i++ )
	{
		allgap = 0;
		for( j=0; j<on; j++ )
		{
			tmpc = oseq[j][i];
			if( tmpc != '-' && tmpc != '=' ) break;
		}
		if( j == on ) 
			allgap = 1;

		if( allgap )
		{
			eqseq[i] = '=';
		}
		else
		{
			eqseq[i] = 'o';
		}
	}

//	for( j=0; j<1; j++ ) reporterr( "\no                = %s\n", oseq[j] );
//	reporterr( "\ne                = %s\n", eqseq );
//	for( j=0; j<1; j++ ) reporterr( "a                = %s\n", aseq[j] );

	if( deletelist )
	{
		for( j=0; j<an; j++ )
		{
			ndel = 0;
			for( i=0,q=0; i<len; i++ )
			{
				tmpc = aseq[j][i];
				if( tmpc != '-' && tmpc != '=' ) 
				{
					if( eqseq[i] == '=' )
					{
//						reporterr( "deleting %d-%d, %c\n", j, i, aseq[j][i] );
						deletelist[j] = realloc( deletelist[j], sizeof( int ) * (ndel+2) );
						deletelist[j][ndel] = q;
						ndel++;
					}
					q++;
				}
			}
			deletelist[j][ndel] = -1;
		}
	}
	for( i=0,p=0; i<len; i++ )
	{

//		if( oseq[0][i] != '=' )
//		reporterr( "i=%d, p=%d, q=%d, originally, %c\n", i, p, q, originallygapped[p]);
//		if( eqseq[i] != '=' && originallygapped[p] != '-' ) // dame!!
		if( eqseq[i] != '=' )
		{
//			reporterr( "COPY! p=%d\n", p );
			if( p != i )
			{
				for( j=0; j<on; j++ ) oseq[j][p] = oseq[j][i];
				for( j=0; j<an; j++ ) aseq[j][p] = aseq[j][i];
			}
			p++;
		}
	}
//		reporterr( "deletemap        = %s\n", deletemap );
//		reporterr( "eqseq            = %s\n", eqseq );
//		reporterr( "originallygapped = %s\n", originallygapped );
	for( j=0; j<on; j++ ) oseq[j][p] = 0;
	for( j=0; j<an; j++ ) aseq[j][p] = 0;

	free( eqseq );

//	for( j=0; j<on; j++ ) reporterr( "\no=%s\n", oseq[j] );
//	for( j=0; j<an; j++ ) reporterr( "a=%s\n", aseq[j] );

	return( i-p );
}

int deletenewinsertions_whole( int on, int an, char **oseq, char **aseq, int **deletelist )
{
	int i, j, p, q, allgap, ndel;
	int len = strlen( oseq[0] );
	char *eqseq, tmpc;

//	reporterr( "In deletenewinsertions_whole\n" );
//	for( j=0; j<on; j++ ) reporterr( "\no=%s\n", oseq[j] );
//	for( j=0; j<an; j++ ) reporterr( "a=%s\n", aseq[j] );

	eqseq = calloc( len+1, sizeof( char ) );
	for( i=0,p=0; i<len; i++ )
	{
		allgap = 0;
		for( j=0; j<on; j++ )
		{
			tmpc = oseq[j][i];
			if( tmpc != '-' ) break;
		}
		if( j == on ) 
			allgap = 1;

		if( allgap )
		{
			eqseq[i] = '=';
		}
		else
		{
			eqseq[i] = 'o';
		}
	}

//	for( j=0; j<1; j++ ) reporterr( "\no                = %s\n", oseq[j] );
//	reporterr( "\ne                = %s\n", eqseq );
//	for( j=0; j<1; j++ ) reporterr( "a                = %s\n", aseq[j] );

	if( deletelist )
	{
		for( j=0; j<an; j++ )
		{
			ndel = 0;
			for( i=0,q=0; i<len; i++ )
			{
				tmpc = aseq[j][i];
				if( tmpc != '-' ) 
				{
					if( eqseq[i] == '=' )
					{
//						reporterr( "deleting %d-%d, %c\n", j, i, aseq[j][i] );
						deletelist[j] = realloc( deletelist[j], sizeof( int ) * (ndel+2) );
						deletelist[j][ndel] = q;
						ndel++;
					}
					q++;
				}
			}
			deletelist[j][ndel] = -1;
		}
	}
	for( i=0,p=0; i<len; i++ )
	{

//		if( oseq[0][i] != '=' )
//		reporterr( "i=%d, p=%d, q=%d, originally, %c\n", i, p, q, originallygapped[p]);
//		if( eqseq[i] != '=' && originallygapped[p] != '-' ) // dame!!
		if( eqseq[i] != '=' )
		{
//			reporterr( "COPY! p=%d\n", p );
			if( p != i )
			{
				for( j=0; j<on; j++ ) oseq[j][p] = oseq[j][i];
				for( j=0; j<an; j++ ) aseq[j][p] = aseq[j][i];
			}
			p++;
		}
	}
//		reporterr( "deletemap        = %s\n", deletemap );
//		reporterr( "eqseq            = %s\n", eqseq );
//		reporterr( "originallygapped = %s\n", originallygapped );
	for( j=0; j<on; j++ ) oseq[j][p] = 0;
	for( j=0; j<an; j++ ) aseq[j][p] = 0;

	free( eqseq );

//	for( j=0; j<on; j++ ) reporterr( "\no=%s\n", oseq[j] );
//	for( j=0; j<an; j++ ) reporterr( "a=%s\n", aseq[j] );
	return( i-p );

}


#if 0
int maskoriginalgaps( char *repseq, char *originallygapped )
{
	int i, p;
	int len = strlen( repseq );
//	reporterr( "repseq = %s\n", repseq );
	for( i=0,p=0; i<len; i++ )
	{
		if( repseq[i] == '=' )
		{
			if( originallygapped[p] == '-' )
			{
				repseq[i] = '-';
				p++;
			}
		}
		else
		{
			p++;
		}
	}
	reporterr( "repseq = %s\n", repseq );
exit( 1 );
}

void restoregaponlysites( char *originallygapped, int n1, int n2, char **s1, char **s2, int rep )
{
	int i, j, p;
	char *tmpnew;
	int len;
	reporterr( "originallygapped = %s\n", originallygapped );
	reporterr( "s1[0]            = %s\n", s1[0] );
	reporterr( "s1[rep]          = %s\n", s1[rep] );
	reporterr( "s2[0]            = %s\n", s2[0] );
exit( 1 );

	tmpnew = calloc( strlen( originallygapped )+1, sizeof( char ) );
	len = strlen( s1[0] );

	for( i=0,p=0; i<len; i++ )
	{
		reporterr( "i=%d, p=%d, s[]=%c, o[]=%c\n", i, p, s1[0][i], originallygapped[p] );
		if( originallygapped[p] == 'o' )
		{
			tmpnew[p] = s1[0][i];
			p++;
		}
		while( originallygapped[p] == '-' )
		{
			tmpnew[p] = '-';
			p++;
		}
	}
	reporterr( "s1[0]            = %s\n", s1[0] );
	reporterr( "tmpnew           = %s\n", tmpnew );
	
}

#endif


int recordoriginalgaps( char *originallygapped, int n, char **s )
{
	int i, j;
	int len = strlen( s[0] );
	int v = 0;
	for( i=0; i<len; i++ )
	{
		for( j=0; j<n; j++ ) if( s[j][i] != '-' ) break;

		if( j == n ) 
			originallygapped[i] = '-';
		else
			originallygapped[i] = 'o';
	}
	originallygapped[i] = 0;
	return( v );
}

void restoreoriginalgaps( int n, char **seq, char *originalgaps )
{
	int i, j, p;
	int lenf = strlen( originalgaps );
	char *tmpseq = calloc( lenf+1, sizeof( char ) );

	for( i=0; i<n; i++ )
	{
		for( j=0,p=0; j<lenf; j++ )
		{
			if( originalgaps[j] == '-' )
				tmpseq[j] = '-';
			else
				tmpseq[j] = seq[i][p++];
		}
		strcpy( seq[i], tmpseq );
	}
	free( tmpseq );
}

void reconstructdeletemap( int nadd, char **addbk, int **deletelist, char **realn, FILE *fp, char **name )
{
	int i, j, p, len;
	char *gapped, *nameptr, *tmpptr;

	for( i=0; i<nadd; i++ )
	{
		len = strlen( addbk[i] );
		gapped = calloc( len+1, sizeof( char ) );
//		for( j=0; j<len; j++ ) gapped[j] = 'o'; // iranai
//		gapped[len] = 0; // iranai

		nameptr = name[i] + 1;
		if( outnumber )
			nameptr = strstr( nameptr, "_numo_e" ) + 8;

		if( (tmpptr=strstr( nameptr, "_oe_" )) ) nameptr = tmpptr + 4; // = -> _ no tame

		fprintf( fp, ">%s\n", nameptr );
		fprintf( fp, "# letter, position in the original sequence, position in the reference alignment\n" );

//		reporterr( "addbk[%d] = %s\n", i, addbk[i] );
		for( j=0; (p=deletelist[i][j])!=-1; j++ )
		{
//			reporterr( "deleting %d, %c\n", p, addbk[i][p] );
			gapped[p] = '-';
		}

//		reporterr( "addbk  = %s\n", addbk[i] );
//		reporterr( "gapped = %s\n", gapped );

		for( j=0,p=0; j<len; j++ )
		{
			while( realn[i][p] == '-' )
				p++;

			if( gapped[j] == '-' )
			{
				fprintf( fp, "%c, %d, -\n", addbk[i][j], j+1 ); // 1origin
			}
			else
			{
				fprintf( fp, "%c, %d, %d\n", addbk[i][j], j+1, p+1 ); // 1origin
				p++;
			}
		}
		free( gapped );
	}
}
