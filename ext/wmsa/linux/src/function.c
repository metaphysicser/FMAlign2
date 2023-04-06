#include "msa.h"
#include "io.h"
#include <spawn.h>
#include <sys/wait.h>
#include <stdlib.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>

// command function

extern char **environ;

void inttostr(char buffer[], int x)
{
    sprintf(buffer, "%d", x);
}

void doubletostr(char buffer[], double x)
{
    sprintf(buffer, "%f", x);
}

int system_spawn(char *prog, char **argv, posix_spawn_file_actions_t *action)
{
#if PRINTARG
    int i;
    for(i = 0; i < 30; ++ i) 
        if(argv[i]) puts(argv[i]); 
        else break;
#endif
    pid_t pid;
    int status = posix_spawn(&pid, prog, action, NULL, argv, NULL);
    // stop until the subprogram ended
    if(! status)
    {
        if (waitpid(pid, &status, 0) == -1)
        {
            reporterr("waitpid");
        }
    }
    else
    {
        reporterr("Error on posix_spawn: %s, calling %s\n", strerror(status), prog);
        exit(1);
    }
    return status;
}

// You can use another cluster algorithm by modifing this function. 
int clustercommand(char *input, char *progfolder)
{
    char *argv[100] = {NULL};
    char *conststr[] = {"-i", "-o", "-c", "-T", "-M", ">", "cd-hit"};
    char buffer[10][10], tmpfilename[maxf], program[maxf];
    int argc = 0, buffers = 0, exitval;
    posix_spawn_file_actions_t action;
    posix_spawn_file_actions_init(&action);
    if(seq_type == 1) // DNA
    {
        sprintf(program, "%s/cd-hit/cd-hit-est", progfolder);
        if(access(program, X_OK)) 
        {
            reporterr("The program 'cd-hit-est' in cd-hit folder is not exist. Please run make all in your program folder or reinstall this program.\n");
            exit(1);
        }
        // sprintf(cmdstr, "-i %s -o %stmp", input, tmpdir);        
    }
    else
    {
        sprintf(program, "%s/cd-hit/cd-hit", progfolder);
        if(access(program, X_OK)) 
        {
            reporterr("The program 'cd-hit' in cd-hit folder is not exist. Please run make all in your program folder or reinstall this program.\n");
            exit(1);
        }
        // sprintf(cmdstr, "-i %s -o %stmp", input, tmpdir);
    }
    // -i input -o tmpdir/tmp
    argv[argc ++] = conststr[6];
    argv[argc ++] = conststr[0];
    argv[argc ++] = input;
    argv[argc ++] = conststr[1];
    sprintf(tmpfilename, "%stmp", tmpdir);
    argv[argc ++] = tmpfilename;
    if(cdhitsim != NOTKNOWNDOUBLE) 
    { 
        // sprintf(cmdstr, "%s -c %f", cmdstr, cdhitsim);
        argv[argc ++] = conststr[2];
        doubletostr(buffer[buffers], cdhitsim);
        argv[argc ++] = buffer[buffers ++];
    }
    // sprintf(cmdstr, "%s -T %d", cmdstr, threads);
    argv[argc ++] = conststr[3];
    inttostr(buffer[buffers], threads);
    argv[argc ++] = buffer[buffers ++];
    // sprintf(cmdstr, "%s -M %d", cmdstr, maxmemory);
    argv[argc ++] = conststr[4];
    inttostr(buffer[buffers], maxmemory);
    argv[argc ++] = buffer[buffers ++];
    if(! printdebug) 
    {
        // sprintf(cmdstr, "%s > /dev/null", cmdstr);
        posix_spawn_file_actions_addopen(&action, STDOUT_FILENO, "/dev/null", O_WRONLY|O_APPEND, 0);
    }
    argv[argc] = NULL;
    exitval = system_spawn(program, argv, &action);
    posix_spawn_file_actions_destroy(&action);
    return exitval;
} 

int maketmpfolder()
{
    sprintf(cmdstr, "mkdir -p %s", tmpdir);
    char *argv[] = {"sh", "-c", cmdstr, NULL};
    reporterr("Info: The folder %s is in your now directory.\n", tmpdir);
    return system_spawn("/bin/sh", argv, NULL);
}

int movefile()
{
    char *argv[] = {"sh", "-c", cmdstr, NULL};
    reporterr("Info: only one cluster in cluster alignment, moved the file into result.\n");
    return system_spawn("/bin/sh", argv, NULL);
}

void cleantmpfile()
{
    sprintf(cmdstr, "rm -rf %s tmp", tmpdir);
    char *argv[] = {"sh", "-c", cmdstr, NULL};
    int val = system_spawn("/bin/sh", argv, NULL);
    if(val) 
    {
        reporterr("Error: tmpdir %s folder cannot remove. You need to remove this folder manually.\n", tmpdir);
    }
}

int calcuate_cluster_sequences()
{
    sprintf(cmdstr, "tmp");
    // reporterr("%s\n", input);
    FILE *cluster = fopen(cmdstr, "r");
    int res;
    int fscanfval = fscanf(cluster, "%d", &res);
    if(fscanfval != 1)
    {
        reporterr("Error on processing file: tmp\n");
        exit(1);
    }
    fclose(cluster);
    return res;
}

// check the existance of BLASTalign, change it if you want use another BLASTalign algorithm
void checkBLASTalign(char *progfolder, char *progname)
{
    sprintf(cmdstr2, "%s/%s", progfolder, progname); // Only in Linux
    if(access(cmdstr2, X_OK))
    {
        reporterr("The program 'staralign' in mafft folder is not exist. Please run make all in your program folder.\n");
        exit(1);
    }
}

// You can use another BLAST alignment algorithm by modifing this function.
int BLASTaligncommand(char *center, char *common, char *progfolder, char *progname)
{
    char buffer[30][10], res[maxf], program[maxf];
    int argc = 0, buffers = 0, exitval;
    char *argv[100] = {NULL};
    char *conststr[] = {"staralign", "-i", "-c", "-f", "-V", "-D", "-P", "-z", "-w", "-B", "-T", "-s", "-b"};
    posix_spawn_file_actions_t action;
    posix_spawn_file_actions_init(&action);
    argv[argc ++] = conststr[0];
    // sprintf(cmdstr, "-i %s -c %s", common, center);
    argv[argc ++] = conststr[1];
    argv[argc ++] = common;
    argv[argc ++] = conststr[2];
    argv[argc ++] = center;
    if(ppenalty != NOTKNOWNINT) 
    {
        // sprintf(cmdstr, "%s -f %d", cmdstr, ppenalty);
        argv[argc ++] = conststr[3];
        inttostr(buffer[buffers], ppenalty);
        argv[argc ++] = buffer[buffers ++];
    }
    if(ppenalty_dist != NOTKNOWNINT) 
    {
        // sprintf(cmdstr, "%s -V %d", cmdstr, ppenalty_dist);
        argv[argc ++] = conststr[4];
        inttostr(buffer[buffers], ppenalty);
        argv[argc ++] = buffer[buffers ++];
    }
    if(seq_type == 1) 
    {
        // sprintf(cmdstr, "%s -D", cmdstr);
        argv[argc ++] = conststr[5];
    }
    else if(seq_type == 2) 
    {
        // sprintf(cmdstr, "%s -P", cmdstr);
        argv[argc ++] = conststr[6];
    }
    if(fftthreshold != NOTKNOWNINT) 
    {
        // sprintf(cmdstr, "%s -z %d", cmdstr, fftthreshold);
        argv[argc ++] = conststr[7];
        inttostr(buffer[buffers], fftthreshold);
        argv[argc ++] = buffer[buffers ++];
    }
    if(fftWinSize != NOTKNOWNINT) 
    {
        // sprintf(cmdstr, "%s -w %d", cmdstr, fftWinSize);
        argv[argc ++] = conststr[8];
        inttostr(buffer[buffers], fftWinSize);
        argv[argc ++] = buffer[buffers ++];
    }
    if(alignband != NOTKNOWNINT) 
    {
        // sprintf(cmdstr, "%s -B %d", cmdstr, alignband);
        argv[argc ++] = conststr[9];
        inttostr(buffer[buffers], alignband);
        argv[argc ++] = buffer[buffers ++];
    }
    if(threads != NOTKNOWNINT) 
    {
        // sprintf(cmdstr, "%s -T %d", cmdstr, threads);
        argv[argc ++] = conststr[10];
        inttostr(buffer[buffers], threads);
        argv[argc ++] = buffer[buffers ++];
    }
    // sprintf(cmdstr, "%s -s %d", cmdstr, nmax_shift);
    argv[argc ++] = conststr[11];
    inttostr(buffer[buffers], nmax_shift);
    argv[argc ++] = buffer[buffers ++];
    // sprintf(cmdstr, "%s -s %d", cmdstr, BLOSUM);
    if(BLOSUM != NOTKNOWNINT && seq_type == 2)
    {
        argv[argc ++] = conststr[12];
        inttostr(buffer[buffers], BLOSUM);
        argv[argc ++] = buffer[buffers ++];
    }
    // sprintf(cmdstr, "%s > %s.res", cmdstr, common);
    sprintf(res, "%s.res", common);
    posix_spawn_file_actions_addopen(&action, STDOUT_FILENO, res, O_CREAT|O_WRONLY|O_TRUNC, S_IRUSR|S_IWUSR|S_IRGRP|S_IROTH);
    if(! printdebug) 
    {
        // sprintf(cmdstr, "%s 2>/dev/null", cmdstr);
        posix_spawn_file_actions_addopen(&action, STDERR_FILENO, "/dev/null", O_WRONLY|O_APPEND, 0);
    }
    argv[argc] = NULL;
    sprintf(program, "%s/%s", progfolder, progname);
    exitval = system_spawn(program, argv, &action);
    posix_spawn_file_actions_destroy(&action);
    return exitval;
}

// check the existance of profile align, change it if you want use another BLASTalign algorithm
void checkprofilealign(char *progfolder, char *progname)
{
    sprintf(cmdstr2, "%s/%s", progfolder, progname); // Only in Linux
    if(access(cmdstr2, X_OK))
    {
        reporterr("The program 'profilealign' in mafft folder is not exist. Please run make all in your program folder.\n");
        exit(1);
    }
}

// You can use another profile alignment algorithm by modifing this function.
int profilealigncommand(char *filelist, char *center, int alignmode, char *progfolder, char *progname)
{
    char buffer[30][10], res[maxf], program[maxf];
    int argc = 0, buffers = 0, exitval;
    char *argv[100] = {NULL};
    char *conststr[] = {"profilealign", "-i", "-p", "-f", "-V", "-D", "-P", "-z", "-w", "-B", "-T", "-A", "-F", "-b"};
    posix_spawn_file_actions_t action;
    posix_spawn_file_actions_init(&action);

    // sprintf(cmdstr, "-i %s -p %s", filelist, center);
    argv[argc ++] = conststr[0];
    argv[argc ++] = conststr[1];
    argv[argc ++] = filelist;
    argv[argc ++] = conststr[2];
    argv[argc ++] = center;
    if(ppenalty != NOTKNOWNINT) 
    {
        // sprintf(cmdstr, "%s -f %d", cmdstr, ppenalty);
        argv[argc ++] = conststr[3];
        inttostr(buffer[buffers], ppenalty);
        argv[argc ++] = buffer[buffers ++];
    }
    if(ppenalty_dist != NOTKNOWNINT) 
    {
        // sprintf(cmdstr, "%s -V %d", cmdstr, ppenalty_dist);
        argv[argc ++] = conststr[4];
        inttostr(buffer[buffers], ppenalty_dist);
        argv[argc ++] = buffer[buffers ++];
    }
    if(seq_type == 1) 
    {
        // sprintf(cmdstr, "%s -D", cmdstr);
        argv[argc ++] = conststr[5];
    }
    else if(seq_type == 2) 
    {
        // sprintf(cmdstr, "%s -P", cmdstr);
        argv[argc ++] = conststr[6];
    }
    if(fftthreshold != NOTKNOWNINT) 
    {
        // sprintf(cmdstr, "%s -z %d", cmdstr, fftthreshold);
        argv[argc ++] = conststr[7];
        inttostr(buffer[buffers], fftthreshold);
        argv[argc ++] = buffer[buffers ++];
    }
    if(fftWinSize != NOTKNOWNINT) 
    {
        // sprintf(cmdstr, "%s -w %d", cmdstr, fftWinSize);
        argv[argc ++] = conststr[8];
        inttostr(buffer[buffers], fftWinSize);
        argv[argc ++] = buffer[buffers ++];
    }
    if(alignband != NOTKNOWNINT) 
    {
        // sprintf(cmdstr, "%s -B %d", cmdstr, alignband);
        argv[argc ++] = conststr[9];
        inttostr(buffer[buffers], alignband);
        argv[argc ++] = buffer[buffers ++];
    }
    if(threads != NOTKNOWNINT) 
    {
        // sprintf(cmdstr, "%s -T %d", cmdstr, threads);
        argv[argc ++] = conststr[10];
        inttostr(buffer[buffers], threads);
        argv[argc ++] = buffer[buffers ++];
    }
    if(alignmode == 0) 
    {
        sprintf(cmdstr, "%s -A ", cmdstr);
        argv[argc ++] = conststr[11];
    }
    if(alignmode == 1) 
    {
        sprintf(cmdstr, "%s -F ", cmdstr);
        argv[argc ++] = conststr[12];
    }
    if(BLOSUM != NOTKNOWNINT && seq_type == 2)
    {
        argv[argc ++] = conststr[13];
        inttostr(buffer[buffers], BLOSUM);
        argv[argc ++] = buffer[buffers ++];
    }
    // sprintf(cmdstr, "%s > %s", cmdstr, outputfile);
    posix_spawn_file_actions_addopen(&action, STDOUT_FILENO, outputfile, O_CREAT|O_WRONLY|O_TRUNC, S_IRUSR|S_IWUSR|S_IRGRP|S_IROTH);
    if(! printdebug) 
    {
        // sprintf(cmdstr, "%s 2>/dev/null", cmdstr);
        posix_spawn_file_actions_addopen(&action, STDERR_FILENO, "/dev/null", O_WRONLY|O_APPEND, 0);
    }
    sprintf(program, "%s/%s", progfolder, progname);
    exitval = system_spawn(program, argv, &action);
    posix_spawn_file_actions_destroy(&action);
    return exitval;
}

static void commongappickpair( char *r1, char *r2, char *i1, char *i2 )
{
	while( *i1 )
	{
		if( *i1 == '-' && *i2 == '-' ) 
		{
			i1++;
			i2++;
		}
		else
		{
			*r1++ = *i1++;
			*r2++ = *i2++;
		}
	}
	*r1 = 0;
	*r2 = 0;
}


double naivepairscore11( char *seq1, char *seq2, int penal, int seq_type )
{
	double  vali;
	int len = strlen( seq1 ), i, j;
	char *s1, *s2, *p1, *p2;
	s1 = AllocateCharVec(len + 10);
	s2 = AllocateCharVec(len + 10);
    vali = 0.0;
    commongappickpair( s1, s2, seq1, seq2 );
    p1 = s1; p2 = s2;
    while( *p1 )
    {
        if( *p1 == '-' )
        {
            vali += (double)penal;
            while( *p1 == '-' )
            {
                p1++;
                p2++;
            }
            continue;
        }
        if( *p2 == '-' )
        {
            vali += (double)penal;
            while( *p2 == '-' )
            {
                p1++;
                p2++;
            }
            continue;
        }
        if(seq_type == 2) 
        {
            for(i = 0; i < 20; ++ i) if(orderprotein[i] == *p1) break;
            for(j = 0; j < 20; ++ j) if(orderprotein[j] == *p2) break;
            ++ p1, ++ p2;
            vali += (double)BLOSUM62[i][j];
        }
        if(seq_type == 1)
        {
            for(i = 0; i < 4; ++ i) if(orderDNA[i] == *p1 || (i == 2 && *p1 == 'u')) break;
            for(j = 0; j < 4; ++ j) if(orderDNA[j] == *p1 || (j == 2 && *p2 == 'u')) break;
            ++ p1, ++ p2;
            vali += (double)trans[i][j];
        }
    }
	if(s1) free( s1 );
	if(s2) free( s2 );
	return( vali );
}


char *get_exe_path(char *buf, int count)
{
    int i;
    int rslt = readlink("/proc/self/exe", buf, count - 1);
    if (rslt < 0 || (rslt >= count - 1))
    {
        return NULL;
    }
    buf[rslt] = 0;
    for (i = rslt; i >= 0; i--)
    {
        if (buf[i] == '/' || buf[i] == '\\')  // '/' in Linux, '\\' in Windows
        {
            buf[i + 1] = 0;
            break;
        }
    }
    return buf;
}
