#include "io.h"
#include "msa.h"

int myatoi( char *in, char before )
{
	if( in == NULL )
	{
		fprintf( stderr, "No argument on -%c.\nPlease check your command, or type \'wmsa -H\' for more help.\n", before );
		exit( 1 );
	}
	return( atoi( in ) );
}

double myatof( char *in, char before )
{
	if( in == NULL )
	{
		fprintf( stderr, "No argument on -%c.\nPlease check your command, or type \'wmsa -H\' for more help.\n", before );
		exit( 1 );
	}
	return( atof( in ) );
}

void reporterr( const char *str, ... )
{
	va_list args;
	va_start( args, str );
	vfprintf( stderr, str, args );
	va_end( args );
	return;
}

int myfgets(char s[], int l, FILE *fp)
{
        int c = 0, i = 0;

		if( feof( fp ) ) return( 1 );

		for( i=0; i<l && ( c=getc( fp ) ) != '\n'; i++ ) 
        	*s++ = c;
        *s = '\0' ;
		if( c != '\n' ) 
			while( getc(fp) != '\n' )
				;
		return( 0 );
}

int countKUorWA( FILE *fp )
{
	int value;
	int c, b;

	value = 0;
	b = '\n';
	while( ( c = getc( fp ) ) != EOF )
	{
		if( b == '\n' && ( c == '>' ) )
			value++;
		b = c;
	}
	rewind( fp );
	return( value );
}

void searchKUorWA( FILE *fp )
{
	int c, b;
	b = '\n';
	while( !( ( ( c = getc( fp ) ) == '>' || c == EOF ) && b == '\n' ) )
		b = c;
	ungetc( c, fp );
}

int countATGC( char *s, int *total )
{
	int nATGC;
	int nChar;
	char c;
	nATGC = nChar = 0;

	if( *s == 0 ) 
	{
		*total = 0;
		return( 0 );
	}

	do
	{
		c = tolower( *s );
		if( isalpha( c ) )
		{
			nChar++;
			if( c == 'a' || c == 't' || c == 'g' || c == 'c' || c == 'u' || c == 'n' )
				nATGC++;
		}
	}
	while( *++s );

	*total = nChar;
	return( nATGC );
}

char *load1SeqWithoutName_realloc( FILE *fpp )
{
	int c, b;
	char *cbuf;
	int size = maxn;
	char *val;

	val = malloc( (size+1) * sizeof( char ) );
	cbuf = val;

	b = '\n';
	while( ( c = getc( fpp ) ) != EOF &&           
          !( ( c == '>' || c == EOF ) && b == '\n' ) )
	{
		*cbuf++ = (char)c;
		if( cbuf - val == size )
		{
			size += N;
			fprintf( stderr, "reallocating...\n" );
			val = (char *)realloc( val, (size+1) * sizeof( char ) );
			if( !val )
			{
				fprintf( stderr, "Allocation error in load1SeqWithoutName_realloc \n" );
				exit( 1 );
			}
			fprintf( stderr, "done.\n" );
			cbuf = val + size - maxn;
		}
		b = c;
	}
	ungetc( c, fpp );
	*cbuf = 0;

	return( val );
}

int getnumlenandtype(FILE *fp, int seq_type)
{
	int total;
	int nsite = 0;
	int atgcnum;
	int i, tmp;
	char *tmpseq;
	char *tmpname;
	double atgcfreq;

#if mingw
	setmode( fileno( fp ), O_BINARY );
	setmode( fileno( stdout ), O_BINARY );
#endif
	seq_maxlen = 0;
	tmpname = AllocateCharVec( N );
	seq_num = countKUorWA( fp );
	searchKUorWA( fp );
	atgcnum = 0;
	total = 0;
	for( i=0; i < seq_num; ++ i )
	{
		myfgets( tmpname, N - 1, fp );
		tmpseq = load1SeqWithoutName_realloc( fp );
		tmp = strlen( tmpseq );
		if(tmp > seq_maxlen) seq_maxlen = tmp;
		atgcnum += countATGC( tmpseq, &nsite );
		total += nsite;
		free( tmpseq );
	}

	atgcfreq = (double)atgcnum / total;
	if(seq_type == 0)
	{
		if( atgcfreq > 0.75 ) seq_type = 1;
		else seq_type = 2;
	}
	free( tmpname );
	return seq_type;
}

void readData( FILE *fp, char **name, int *nlen, char **seq, int njob )
{
	int i; 
	static char *tmpseq = NULL;

	rewind( fp );
	searchKUorWA( fp );

	for( i=0; i<njob; i++ )
	{
		name[i][0] = '='; getc( fp ); 
#if 0
		fgets( name[i]+1, B-2, fp ); 
		j = strlen( name[i] );
		if( name[i][j-1] != '\n' )
			ErrorExit( "Too long name\n" );
		name[i][j-1] = 0;
#else
		myfgets( name[i]+1, B-2, fp ); 
#endif
		tmpseq = load1SeqWithoutName_realloc( fp );
		strcpy( seq[i], tmpseq );
		free( tmpseq );
		nlen[i] = strlen( seq[i] );
	}
}

void gapfilter_oneseq( char *aseq, char *seq )
{
	for( ; *seq != 0; seq++ )
	{
		if( *seq != '-' )
			*aseq++ = *seq;
	}
	*aseq = 0;
}

void writeData( FILE *fp, int locnjob, char **name, int *nlen, char **aseq )
{
	int i, j;
	int nalen;

	for( i=0; i<locnjob; i++ )
	{
		nalen = strlen( aseq[i] );
		fprintf( fp, ">%s\n", name[i]+1 );
		fprintf(fp, "%s\n", aseq[i]);
	}
}
