// command function

int clustercommand(char *input, char *progfolder); // cluster program processor
int maketmpfolder(); // make temp folder, exit value is not 0 if failed on making this folder
void cleantmpfile(); // clear all the temp file in the program
int calcuate_cluster_sequences(); // calcuate the number of cluster sequences
void checkBLASTalign(char *progfolder, char *programname); // check the existance BLASTalign program
int BLASTaligncommand(char *center, char *common, char *progfolder, char *progname); // staralign command processor
void checkprofilealign(char *progfolder, char *programname); // check the existance profile align program
int profilealigncommand(char *filelist, char *center, int alignmode, char *progfolder, char *progname); // profilealign command processor
double naivepairscore11( char *seq1, char *seq2, int penal, int seq_type); // SP Scores (one by one)
char *get_exe_path(char *buf, int count); // find the program path
int movefile(); // move file, the command is written on the cmdstr
