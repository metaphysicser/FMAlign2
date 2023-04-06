#include "mltaln.h"
#include "threadpool.h"
#include <time.h>

#define REPORTCOSTS 1
#define SHOWPROFILEVERSION reporterr( "%s (%s) Version " VERSION "\n%d thread(s)\n\n", progName( argv[0] ), (dorp=='d')?"nuc":((nblosum==-2)?"text":"aa"), nthread )

#define DEBUG 0
#define TESTCONST 0

#define ITERATIVECYCLE 2

#define END_OF_VEC -1

static int tuplesize, nunknown, exitval, aligncases;

#define PLENFACA 0.01
#define PLENFACB 10000
#define PLENFACC 10000
#define PLENFACD 0.1
#define D6LENFACA 0.01
#define D6LENFACB 2500
#define D6LENFACC 2500
#define D6LENFACD 0.1
#define D10LENFACA 0.01
#define D10LENFACB 1000000
#define D10LENFACC 1000000
#define D10LENFACD 0.0

#ifdef enablemultithread
typedef struct _merge_file_arg
{
	int f1, f2, floor;
	char *f1name, *f2name;
	pthread_mutex_t *deplock;
	pthread_cond_t *cond;
	Treedep *dep;
} merge_file_arg;
#endif

void print_help()
{
	reporterr("Profile alignment Version %s help:\n", VERSION);
	reporterr("-i: sequences file name, every line has a file name without spaces\n");
	reporterr("-p: center file with FASTA format\n");
	reporterr("-T: use T threads to run this program\n");
	reporterr("-f, -g, -h: ppenalty, ppenalty_ex(not used), poffset(not used)\n");
	reporterr("-Q, -V: penalty_shift_factor(not used), ppenalty_dist\n");
	reporterr("-b: BLOSUM Matrix\n");
	reporterr("-j: use jtt/kimura model, pamN is needed\n");
	reporterr("-m: use tm model, pamN is needed\n");
	reporterr("-D, -P: -D the sequence is DNA, -P the sequence is Protein\n");
	reporterr("-S: Calcuate SP Scores after alignment\n");
	reporterr("-z, -w: FFT align arguments: fftthreshold, fftWinSize\n");
	reporterr("-B: Kband in calcuating DP-matrix during the alignment\n");
	reporterr("-W: tuplesize on UPGMA (6 or 10, 10 is only for DNA/RNA sequences)\n");
	reporterr("-X: Use mix method to calcuate the UPGMA cluster, the sub effect is needed, it must be in (0, 1)\n");
	reporterr("-E: Use average method to calcuate the UPGMA cluster\n");
	reporterr("-q: Use minimum method to calcuate the UPGMA cluster\n");
	reporterr("-A: Use Aalign to align sequences\n");
	reporterr("-F: Use FFT align to align sequences\n");
	reporterr("-v: show program version and exit\n");
	reporterr("-H, -?: Print help message and exit\n");
}

void print_version()
{
	reporterr("profilealign %s\n", VERSION);
}

void arguments( int argc, char *argv[] )
{
	int c;
	alg = 'A';
	nthread = 1;
	outnumber = 0;
	nevermemsave = 0;
	inputfile = NULL;
	nblosum = 62;
	treemethod = 'X';
	sueff_global = 0.1;
	contin = 0;
	scoremtx = 1;
	kobetsubunkatsu = 0;
	dorp = NOTSPECIFIED;
	ppenalty = -1530;
	penalty_shift_factor = 1000.0;
	poffset = -123;
	kimuraR = NOTSPECIFIED;
	pamN = NOTSPECIFIED;
	fftWinSize = NOTSPECIFIED;
	fftThreshold = NOTSPECIFIED;
	TMorJTT = JTT;
	spscoreout = 0;
	tuplesize = 6;
	aligncases = 1;
	alignband = NOTSPECIFIED;
	
	while( --argc > 0 && (*++argv)[0] == '-' )
	{
		while ( (c = *++argv[0]) )
		{
			switch( c )
			{
				case 'i': // common sequence file per line, must not have space(s)
					inputfile = *++argv;
					reporterr( "file list file = %s\n", inputfile );
					--argc;
					goto nextoption;
				case 'p':
					centerfile = *++ argv;
					reporterr("center sequence file = %s\n", centerfile );
					-- argc;
					goto nextoption;
				case 'f':
					ppenalty = (int)( myatof( *++argv ) * 1000 - 0.5 );
					reporterr(       "ppenalty = %d\n", ppenalty );
					--argc;
					goto nextoption;
				case 'Q':
					penalty_shift_factor = myatof( *++argv );
					--argc;
					goto nextoption;
				case 'g':
					ppenalty_ex = (int)( myatof( *++argv ) * 1000 - 0.5 );
					reporterr(       "ppenalty_ex = %d\n", ppenalty_ex );
					--argc;
					goto nextoption;
				case 'h':
					poffset = (int)( myatof( *++argv ) * 1000 - 0.5 );
//					reporterr(       "poffset = %d\n", poffset );
					--argc;
					goto nextoption;
				case 'k':
					kimuraR = myatoi( *++argv );
					reporterr(       "kimura model, kappa = %d\n", kimuraR );
					--argc;
					goto nextoption;
				case 'b':
					nblosum = myatoi( *++argv );
					reporterr(       "blosum %d / kimura 200 \n", nblosum );
					--argc;
					goto nextoption;
				case 'j':
					pamN = myatoi( *++argv );
					scoremtx = 0;
					TMorJTT = JTT;
					reporterr(       "jtt/kimura %d\n", pamN );
					--argc;
					goto nextoption;
				case 'm':
					pamN = myatoi( *++argv );
					scoremtx = 0;
					TMorJTT = TM;
					reporterr(       "tm %d\n", pamN );
					--argc;
					goto nextoption;
				case 'T':
					nthread = myatoi( *++argv );
					reporterr(       "nthread = %d\n", nthread );
					-- argc; 
					goto nextoption;
				case 'D':
					dorp = 'd';
					break;
				case 'P':
					dorp = 'p';
					break;
				case 'X':
					treemethod = 'X';
					sueff_global = atof( *++argv );
					fprintf( stderr, "sueff_global = %f\n", sueff_global );
					--argc;
					goto nextoption;
  				case 'E':
					treemethod = 'E';
					break;
				case 'q':
					treemethod = 'q'; // minimum
					break;
				case 'S' :
					spscoreout = 1; // 2014/Dec/30, sp score
					break;
				case 'z':
					fftThreshold = myatoi( *++argv );
					--argc; 
					goto nextoption;
				case 'w':
					fftWinSize = myatoi( *++argv );
					--argc;
					goto nextoption;
				case 'W':
					tuplesize = myatoi( *++argv );
					--argc;
					goto nextoption;
				case 'B':
					alignband = myatoi( *++argv );
					-- argc;
					goto nextoption;
				case 'A':
					aligncases = 1; // A__align11
					reporterr("Use Aalign\n");
					break;
				case 'F':
					aligncases = 0; // Falign
					reporterr("Use FFT Align\n");
					break;
				case 'H':
				case '?':
					print_help();
					exit(0);
				case 'v':
					print_version();
					exit(0);
				default:
					reporterr(       "illegal option %c\n", c );
					argc = 0;
					break;
			}
		}
		nextoption:
			;
	}
	if( argc != 0 ) 
	{
		reporterr( "options: Check source file !\n" );
		exit( 1 );
	}
}

void makepointtable_nuc( int *pointt, int *n )
{
	int point;
	register int *p;

	if( *n == -1 )
	{
		*pointt = -1;
		return;
	}

	p = n;
	point  = *n++ *  1024;
	point += *n++ *   256;
	point += *n++ *    64;
	point += *n++ *    16;
	point += *n++ *     4;
	point += *n++;
	*pointt++ = point;

	while( *n != END_OF_VEC )
	{
		point -= *p++ * 1024;
		point *= 4;
		point += *n++;
		*pointt++ = point;
	}
	*pointt = END_OF_VEC;
}

void makepointtable( int *pointt, int *n )
{
	int point;
	register int *p;

	if( *n == -1 )
	{
		*pointt = -1;
		return;
	}

	p = n;
	point  = *n++ *  7776;
	point += *n++ *  1296;
	point += *n++ *   216;
	point += *n++ *    36;
	point += *n++ *     6;
	point += *n++;
	*pointt++ = point;

	while( *n != END_OF_VEC )
	{
		point -= *p++ * 7776;
		point *= 6;
		point += *n++;
		*pointt++ = point;
	}
	*pointt = END_OF_VEC;
}

void makepointtable_nuc_dectet( int *pointt, int *n )
{
	int point;
	register int *p;

	if( *n == -1 )
	{
		*pointt = -1;
		return;
	}

	p = n;
	point  = *n++ *262144;
	point += *n++ * 65536;
	point += *n++ * 16384;
	point += *n++ *  4096;
	point += *n++ *  1024;
	point += *n++ *   256;
	point += *n++ *    64;
	point += *n++ *    16;
	point += *n++ *     4;
	point += *n++;
	*pointt++ = point;

	while( *n != END_OF_VEC )
	{
		point -= *p++ *262144;
		point *= 4;
		point += *n++;
		*pointt++ = point;

	}
	*pointt = END_OF_VEC;
}

void seq_grp_nuc( int *grp, char *seq )
{
	int tmp;
	int *grpbk = grp;
	while( *seq )
	{
		tmp = amino_grp[(int)*seq++];
		if( tmp < 4 )
			*grp++ = tmp;
		else
			nunknown++;
	}
	*grp = END_OF_VEC;
	if( grp - grpbk < tuplesize )
	{
		*grpbk = -1;
	}
}

void seq_grp( int *grp, char *seq )
{
	int tmp;
	int *grpbk = grp;
	while( *seq )
	{
		tmp = amino_grp[(unsigned char)*seq++];
		if( tmp < 6 )
			*grp++ = tmp;
		else
			nunknown++;
	}
	*grp = END_OF_VEC;
	if( grp - grpbk < 6 )
	{
		*grpbk = -1;
	}
}

void makecompositiontable_p( int *table, int *pointt )
{
	int point;

	while( ( point = *pointt++ ) != END_OF_VEC )
		table[point]++;
}

#ifdef enablemultithread
int finished;
void *merge_multithread(void *arg)
{
#define F(X) (((merge_file_arg *)arg) -> X)
	char *f1name = F(f1name), *f2name = F(f2name);
	pthread_mutex_t *deplock = F(deplock);
	pthread_cond_t *cond = F(cond);
	int f1 = F(f1), f2 = F(f2), floor = F(floor);
	Treedep *dep = F(dep);
#undef F
	int f1seq, f2seq, j, len1, len2, processlen, alloclen, fftlog;
	char **seq, **seq2, **name, **name2, *sgap1, *sgap2, *egap1, *egap2;
	int *nlen, *nlen22;
	double *eff, *eff2;
	pthread_mutex_lock(deplock);
	if(dep[floor].child0 != -1) 
		while(! dep[dep[floor].child0].done) 
			pthread_cond_wait(cond, deplock);
	if(dep[floor].child1 != -1) 
		while(! dep[dep[floor].child1].done) 
			pthread_cond_wait(cond, deplock);
	pthread_mutex_unlock(deplock);
	// reporterr("Task %d: %s, %s\n", floor, f1name, f2name);
	FILE *f1fp = fopen(f1name, "r"), *f2fp = fopen(f2name, "r");
	getnumlen_nocommonnjob(f1fp, &f1seq, &len1);
	rewind(f1fp);
	getnumlen_nocommonnjob(f2fp, &f2seq, &len2);
	rewind(f2fp);
	processlen = MAX(len1, len2);
	seq = AllocateCharMtx(f1seq, processlen << 1);
	name = AllocateCharMtx(f1seq, B);
	nlen = AllocateIntVec(f1seq);
	seq2 = AllocateCharMtx(f2seq, processlen << 1);
	name2 = AllocateCharMtx(f2seq, B);
	nlen22 = AllocateIntVec(f2seq);
	readData_pointer2(f1fp, f1seq, name, nlen, seq);
	readData_pointer2(f2fp, f2seq, name2, nlen22, seq2);
	fclose(f1fp);
	fclose(f2fp);
	// Calcuate the effect value: average of all the sequence
	eff = AllocateDoubleVec(f1seq);
	eff2 = AllocateDoubleVec(f2seq);
	if(f1seq > 0) for(j = 0; j < f1seq; ++ j) eff[j] = 1.0 / f1seq;
	else 
	{ 
		reporterr("Error: the sequence file %s has no sequences! It may coursed by the smaller value of fftWinsize, please make it larger. The arugment of fftWinsize is -w.\n", f1name); 
		exit(1); 
		return (void *)-1;
	}
	if(f2seq > 0) for(j = 0; j < f2seq; ++ j) eff2[j] = 1.0 / f2seq;
	else 
	{ 
		reporterr("Error: the sequence file %s has no sequences! It may coursed by the smaller value of fftWinsize, please make it larger. The arugment of fftWinsize is -w.\n", f2name);
		exit(1); 
		return (void *)-1;
	}
	alloclen = len1 + len2 + 10;
	if(aligncases == 1) 
	{
		sgap1 = AllocateCharVec(f1seq + 10);
		sgap2 = AllocateCharVec(f2seq + 10);
		egap1 = AllocateCharVec(f1seq + 10);
		egap2 = AllocateCharVec(f2seq + 10);
		memset(sgap1, 'o', f1seq * sizeof(char));
		memset(sgap2, 'o', f2seq * sizeof(char));
		memset(egap1, 'o', f1seq * sizeof(char));
		memset(egap2, 'o', f2seq * sizeof(char));
		A__align(n_dis_consweight_multi, penalty, penalty_ex, seq, seq2, eff, eff2, f1seq, f2seq, alloclen, sgap1, sgap2, egap1, egap2, 1, 1);
		free(sgap1);
		free(sgap2);
		free(egap1);
		free(egap2);
	}
	else if(aligncases == 0) 
		Falign(NULL, NULL, n_dis_consweight_multi, seq, seq2, eff, eff2, NULL, NULL, f1seq, f2seq, alloclen, &fftlog, NULL, 0, NULL);
	else 
	{
		ErrorExit("ERROR: aligncases is error. Please check your command.\n");
		return (void *)-1;
	}
	f1fp = fopen(f1name, "w");
	writeData_pointer(f1fp, f1seq, name, nlen, seq);
	writeData_pointer(f1fp, f2seq, name2, nlen22, seq2);
	fclose(f1fp);
	f2fp = fopen(f2name, "w");
	fclose(f2fp);

	pthread_mutex_lock(deplock);
	dep[floor].done = 1;
	pthread_cond_broadcast(cond);
	pthread_mutex_unlock(deplock);

	// reporterr("Task %d done, %d, %d\n", floor, f1, f2);
	FreeDoubleVec(eff);
	FreeDoubleVec(eff2);
	FreeCharMtx(seq);
	FreeCharMtx(name);
	FreeIntVec(nlen);
	FreeCharMtx(seq2);
	FreeCharMtx(name2);
	FreeIntVec(nlen22);
	++ finished;
	if(finished % 10 == 0 || finished == njob - 1) reporterr("\r   %d / %d", finished, njob - 1);
	return NULL;
}
#endif

int main(int argc, char **argv)
{
	int *nlen = NULL;
	char **seq = NULL, **name = NULL, *realfilename = NULL, **realfilename2d = NULL, **seq2 = NULL, **name2 = NULL;
	int maxlen, alloclen, i, j, k, fftlog, centerseqs;
	int f1, f2; // sequence place
	int f1seq, f2seq; // sequence numbers
	int f1len, f2len; // max length of two sequences
	int *grpseq = NULL, **pointt = NULL, *table1 = NULL, ***topol = NULL, *nogaplen = NULL, *nlen22;
	char *tmpseq = NULL, *align1 = NULL, *align2 = NULL;
	double *eff = NULL, **mtx = NULL, **nlen2 = NULL, *eff2 = NULL;
	char *sgap1 = NULL, *sgap2 = NULL, *egap1 = NULL, *egap2 = NULL;
	Treedep *dep; 
	double bunbo, lenfac, longer, shorter;
	char b[B];
	FILE *infp = NULL, *cefp = NULL, *f1fp = NULL, *f2fp = NULL;
	arguments(argc, argv);
	if(inputfile && centerfile)
	{
		infp = fopen(inputfile, "rb");
		if(! infp)
		{
			reporterr("Error: cannot open sequences file %s\n", inputfile);
			exit(1);
		}
		cefp = fopen(centerfile, "rb");
		if(! cefp)
		{
			reporterr("Error: cannot open center sequence file %s\n", centerfile);
			exit(1);
		}
	}
	else
	{
		ErrorExit("In profile align, you must specify the center sequence file and the other sequences file. ");
	}
#if !defined(mingw) && !defined(_MSC_VER)
	setstacksize( (unsigned long long)1024 * 1024 * 1024 ); // topolorder() de ookime no stack wo shiyou.
#endif

	/* Part 1: use center sequence file get a UPGMA tree */
	getnumlen(cefp);
	rewind(cefp);
	if(njob == 1) // case (center sequence is only one sequence): the program is not neeeded, print the profile
	{
		reporterr("There is only one sequence on the center file. \n");
		fclose(cefp);
		infp = fopen(inputfile, "rb");
		realfilename = AllocateCharVec(256);
		exitval = fscanf(infp, "%s", realfilename);
		fclose(infp);
		infp = fopen(realfilename, "rb");
		getnumlen(infp);
		rewind(infp);
		name = AllocateCharMtx(njob, nlenmax + 1);
		nlen = AllocateIntVec(njob);
		seq = AllocateCharMtx(njob, nlenmax + 1);
		readData_pointer(infp, name, nlen, seq);
		fclose(infp);

		writeData_pointer(stdout, njob, name, nlen, seq);
		SHOWPROFILEVERSION;
		FreeIntVec(nlen);
		FreeCharMtx(name);
		FreeCharMtx(seq);
		return 0;
	}
	centerseqs = njob;
	
	if( dorp == 'p' && tuplesize != 6 )
		ErrorExit( "tuplesize must be 6 for amino acid sequence\n" );
	if( dorp == 'd' && tuplesize != 6 && tuplesize != 10 )
		reporterr( "tuplesize must be 6 or 10 for DNA/RNA sequence\n" );

	name = AllocateCharMtx(njob, nlenmax + 1);
	nlen = AllocateIntVec(njob);
	seq = AllocateCharMtx(njob, nlenmax + 1);
	readData_pointer(cefp, name, nlen, seq);
#if TESTCONST
	disp = 1;
#endif
	constants(njob, seq);
#if REPORTCOSTS
		time_t starttime, startclock;
		starttime = time(NULL);
		startclock = clock();
#endif

	/* Part 1.2: calcuate the 6-mer or 10-mer distance */
	reporterr( "\nMaking a distance matrix using the center sequences...\n" );
	mtx = AllocateDoubleHalfMtx(njob);
	tmpseq = AllocateCharVec( nlenmax+1 );
	grpseq = AllocateIntVec( nlenmax+1 );
	pointt = AllocateIntMtx( njob, nlenmax+1 );
	nogaplen = AllocateIntVec(njob);
	if( dorp == 'd' ) tsize = (int)pow( 4, tuplesize );
	else              tsize = (int)pow( 6, 6 );
	if( dorp == 'd' && tuplesize == 6 )
	{
		lenfaca = D6LENFACA;
		lenfacb = D6LENFACB;
		lenfacc = D6LENFACC;
		lenfacd = D6LENFACD;
	}
	else if( dorp == 'd' && tuplesize == 10 )
	{
		lenfaca = D10LENFACA;
		lenfacb = D10LENFACB;
		lenfacc = D10LENFACC;
		lenfacd = D10LENFACD;
	}
	else    
	{
		lenfaca = PLENFACA;
		lenfacb = PLENFACB;
		lenfacc = PLENFACC;
		lenfacd = PLENFACD;
	}
	maxl = 0;
	for( i=0; i<njob; i++ ) 
	{
		gappick0( tmpseq, seq[i] );
		nogaplen[i] = strlen( tmpseq );
		if( nogaplen[i] > maxl ) maxl = nogaplen[i];
		if( dorp == 'd' ) /* nuc */
		{
			seq_grp_nuc( grpseq, tmpseq );
			if( tuplesize == 10 )
				makepointtable_nuc_dectet( pointt[i], grpseq );
			else if( tuplesize == 6 )
				makepointtable_nuc( pointt[i], grpseq );
		}
		else              /* amino */
		{
			seq_grp( grpseq, tmpseq );
			makepointtable( pointt[i], grpseq );
		}
	}
	if( nunknown ) reporterr( "\nWarning: There are %d ambiguous characters on centerfile.\n", nunknown );

	for( i=0; i<njob; i++ )
	{
		table1 = (int *)calloc( tsize, sizeof( int ) );
		if( !table1 ) ErrorExit( "Cannot allocate table1\n" );
		if( i % 100 == 0 )
		{
			reporterr( "\r%5d / %d", i+1, njob );
		}
		makecompositiontable_p( table1, pointt[i] );

		for( j=i; j<njob; j++ ) 
		{
			mtx[i][j-i] = (double)commonsextet_p( table1, pointt[j] );
		} 
		free( table1 ); table1 = NULL;
	}
	
	for( i = 0; i < njob - 1; ++ i)
	{
		for( j = i + 1; j < njob; ++ j ) 
		{
			if( nogaplen[i] > nogaplen[j] )
			{
				longer=(double)nogaplen[i];
				shorter=(double)nogaplen[j];
			}
			else
			{
				longer=(double)nogaplen[j];
				shorter=(double)nogaplen[i];
			}
			lenfac = 1.0 / ( shorter / longer * lenfacd + lenfacb / ( longer + lenfacc ) + lenfaca );
			bunbo = MIN( mtx[i][0], mtx[j][0] );
			if( bunbo == 0.0 )
				mtx[i][j-i] = 2.0;
			else
				mtx[i][j-i] = ( 1.0 - mtx[i][j-i] / bunbo ) * lenfac * 2.0;
		}
	}
	FreeIntVec(nogaplen); nogaplen = NULL; 
	FreeIntMtx( pointt ); pointt = NULL;
	free( grpseq ); grpseq = NULL;
	free( tmpseq ); tmpseq = NULL;
	reporterr("\ndone. \n");
	/* Part 1.3: UPGMA tree */
	reporterr("\nUse center sequence to get an UPGMA tree...");
	dep = (Treedep *)malloc(njob * sizeof(Treedep));
	topol = AllocateIntCub(njob, 2, 0);
	nlen2 = AllocateDoubleMtx(njob, 2);
	fixed_musclesupg_double_realloc_nobk_halfmtx_memsave(njob, mtx, topol, nlen2, dep, 1, 0);
	// fixed_musclesupg_double_treeout(njob, mtx, topol, nlen2, name);
	FreeDoubleHalfMtx( mtx, njob ); mtx = NULL;
	// FreeDoubleMtx(nlen2); 
	FreeCharMtx(name); name = NULL;
	FreeCharMtx(seq); seq = NULL;
	FreeIntVec(nlen); nlen = NULL;
	reporterr("\ndone. \n");


	/* Part 2: Use the tree to align profile */
	/* Part 2.1: read file list to get the matrix density */
	realfilename2d = AllocateCharMtx(centerseqs, 200);
	rewind(infp);
	for(i = 0; i < centerseqs; ++ i)
	{
		exitval = fscanf(infp, "%s", realfilename2d[i]);
		if(exitval != 1) ErrorExit("Input file is error, please check!");
		f1fp = fopen(realfilename2d[i], "rb");
		if(! f1fp)
		{
			reporterr("ERROR: The file %s in input file can't open!\n", realfilename2d[i]);
			exit(1);
		}
		fclose(f1fp);
	}
		
	/* Part 2.2: use changed FFT (changed to calcuate density instead of sequence) to align profile data */
	reporterr("Aligning all the profiles...\n");
#ifndef enablemultithread
	// single thread on profile alignment
	for(i = 0; i < centerseqs - 1; ++ i)
	{
		// only need to use topol[i][0][0] && topol[i][1][0]
		f1 = topol[i][0][0], f2 = topol[i][1][0];
#if 0
		reporterr("%d: %d %d %d %d %f %f\n", i, f1, f2, dep[i].child0, dep[i].child1, nlen2[i][0], nlen2[i][1]);
#endif
		// read two sequences where the place is f1 and f2
		f1fp = fopen(realfilename2d[f1], "rb");
		f2fp = fopen(realfilename2d[f2], "rb");
		getnumlen(f1fp);
		rewind(f1fp);
		f1seq = njob;
		maxlen = nlenmax;
		getnumlen(f2fp);
		rewind(f2fp);
		f2seq = njob;
		maxlen = MAX(maxlen, nlenmax);
		nlenmax = maxlen;
#if DEBUG
		reporterr("%d, %d\n", f1seq, f2seq);
#endif
		seq = AllocateCharMtx(f1seq, nlenmax << 1);
		name = AllocateCharMtx(f1seq, B);
		nlen = AllocateIntVec(f1seq);
		seq2 = AllocateCharMtx(f2seq, nlenmax << 1);
		name2 = AllocateCharMtx(f2seq, B);
		nlen22 = AllocateIntVec(f2seq);
		readData_pointer2(f1fp, f1seq, name, nlen, seq);
		readData_pointer2(f2fp, f2seq, name2, nlen22, seq2);
		fclose(f1fp);
		fclose(f2fp);
		// Calcuate the effect value: average of all the sequence
		eff = AllocateDoubleVec(f1seq);
		eff2 = AllocateDoubleVec(f2seq);
		if(f1seq > 0) for(j = 0; j < f1seq; ++ j) eff[j] = 1.0 / f1seq;
		else { reporterr("Error: the sequence file %s has no sequences! It may coursed by the smaller value of fftWinsize, please make it larger. The arugment of fftWinsize is -w.\n", realfilename2d[f1]); exit(1); }
		if(f2seq > 0) for(j = 0; j < f2seq; ++ j) eff2[j] = 1.0 / f2seq;
		else { reporterr("Error: the sequence file %s has no sequences! It may coursed by the smaller value of fftWinsize, please make it larger. The arugment of fftWinsize is -w.\n", realfilename2d[f2]); exit(1); }
		f1len = nlen[0];
		for(j = 1; j < f1seq; ++ j) f1len = MAX(nlen[j], f1len);
		f2len = nlen22[0];
		for(j = 1; j < f2seq; ++ j) f2len = MAX(nlen22[j], f2len);
		alloclen = f1len + f2len + 10;
		if(aligncases == 1) 
		{
			sgap1 = AllocateCharVec(f1seq + 10);
			sgap2 = AllocateCharVec(f2seq + 10);
			egap1 = AllocateCharVec(f1seq + 10);
			egap2 = AllocateCharVec(f2seq + 10);
			memset(sgap1, 'o', f1seq * sizeof(char));
			memset(sgap2, 'o', f2seq * sizeof(char));
			memset(egap1, 'o', f1seq * sizeof(char));
			memset(egap2, 'o', f2seq * sizeof(char));
			A__align(n_dis_consweight_multi, penalty, penalty_ex, seq, seq2, eff, eff2, f1seq, f2seq, alloclen, sgap1, sgap2, egap1, egap2, 1, 1);
			free(sgap1);
			free(sgap2);
			free(egap1);
			free(egap2);
		}
		else if(aligncases == 0) 
			Falign(NULL, NULL, n_dis_consweight_multi, seq, seq2, eff, eff2, NULL, NULL, f1seq, f2seq, alloclen, &fftlog, NULL, 0, NULL);
		else ErrorExit("ERROR: aligncases is error. Please check your command.\n");
		f1fp = fopen(realfilename2d[f1], "w");
		writeData_pointer(f1fp, f1seq, name, nlen, seq);
		writeData_pointer(f1fp, f2seq, name2, nlen22, seq2);
		fclose(f1fp);
		f2fp = fopen(realfilename2d[f2], "w");
		fclose(f2fp);
		FreeDoubleVec(eff);
		FreeDoubleVec(eff2);
		if(i != centerseqs - 2)
		{
			FreeCharMtx(seq);
			FreeCharMtx(name);
			FreeIntVec(nlen);
			FreeCharMtx(seq2);
			FreeCharMtx(name2);
			FreeIntVec(nlen22);
		}
		reporterr("\r   %d / %d", i + 1, centerseqs - 1);
	}
	reporterr("\ndone. \n");
	/* Part 3A: write the answer file */
	writeData_pointer(stdout, f1seq, name, nlen, seq);
	writeData_pointer(stdout, f2seq, name2, nlen22, seq2);
#else
	pthread_mutex_t deplock;
	pthread_cond_t cond;
	for(i = 0; i < centerseqs; ++ i) dep[i].done = 0;
	pthread_cond_init(&cond, NULL);
	pthread_mutex_init(&deplock, NULL);
	threadpool_t tp;
	threadpool_init(&tp, nthread);
	merge_file_arg *_merge_arg_;
	_merge_arg_ = malloc(sizeof(merge_file_arg) * centerseqs);
#if 0
	reporterr("tree: \n");
	for(i = 0; i < centerseqs - 1; ++ i)
		reporterr("floor %d: topol(%d, %d), dep(%d, %d)\n", i, topol[i][0][0], topol[i][1][0], dep[i].child0, dep[i].child1);
#endif
	finished = 0;
	for(i = 0, j; i < centerseqs - 1; ++ i)
	{
		f1 = topol[i][0][0], f2 = topol[i][1][0];
		_merge_arg_[i].f1 = f1;
		_merge_arg_[i].f2 = f2;
		_merge_arg_[i].f1name = realfilename2d[f1];
		_merge_arg_[i].f2name = realfilename2d[f2];
		_merge_arg_[i].dep = dep;
		_merge_arg_[i].floor = i;
		_merge_arg_[i].cond = &cond;
		_merge_arg_[i].deplock = &deplock;
		threadpool_add_task(&tp, merge_multithread, _merge_arg_ + i);
	}
	threadpool_destroy(&tp);
	pthread_mutex_destroy(&deplock);
	pthread_cond_destroy(&cond);
	free(_merge_arg_);
	reporterr("\ndone. \n");
	/* Part 3B: write the answer file */
	f1fp = fopen(realfilename2d[topol[centerseqs - 2][0][0]], "r");
	Filecopy(f1fp, stdout);
	fclose(f1fp);
#endif

	/* Part end: free and show info */
#if REPORTCOSTS
//		use_getrusage();
		reporterr( "\nprofilealign, real = %f min\n", (float)(time(NULL) - starttime)/60.0 );
		reporterr( "profilealign, user = %f min\n", (float)(clock()-startclock)/CLOCKS_PER_SEC/60);
#endif
	if(spscoreout) 
	{
		if(seq) FreeCharMtx(seq);
		if(name) FreeCharMtx(name);
		if(nlen) FreeIntVec(nlen);
		f1fp = fopen(realfilename2d[topol[centerseqs - 2][0][0]], "rb");
		getnumlen(f1fp); rewind(f1fp);
		seq = AllocateCharMtx(njob, nlenmax << 1);
		name = AllocateCharMtx(njob, B + 10);
		nlen = AllocateIntVec(njob);
		readData_pointer(f1fp, name, nlen, seq);
		fclose(f1fp);
		reporterr("SP Scores = %.6f\n", sumofpairsscore(njob, seq));
		FreeIntVec(nlen);
		FreeCharMtx(name);
		FreeCharMtx(seq);
	}
	SHOWPROFILEVERSION;
#ifndef enablemultithread
	reporterr("...but NOT used multi threads in profilealign.\n\n\n");
#endif
	FreeCharMtx(realfilename2d);
	free(dep);
	FreeIntCub(topol);
	return 0;
}
