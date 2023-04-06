CC = gcc
#CC = icc
CFLAGS = -O3
INSTALLFOLDER = /usr/bin
CDHITFOLDER = $(INSTALLFOLDER)/cd-hit
MAFFTFOLDER = $(INSTALLFOLDER)/submafft
# You can edit this folder in order to install the program to another place

# If in MAC OS X, please not use -fopenmp
NOOPENMP = no

ifeq ($(NOOPENMP), yes)
	NOOPENMPFLAG = openmp=no
endif

# Comment out the above line if your compiler 
# does not support TLS (thread-local strage).
ENABLE_MULTITHREAD = -Denablemultithread
THREADS = 1

ifdef ENABLE_MULTITHREAD
LIBS = -lm  -lpthread
else
LIBS = -lm
endif

ifdef ENABLE_ATOMIC
STDF = -std=c11
else
STDF = -std=c99
endif

MYCFLAGS = $(MNO_CYGWIN) $(ENABLE_MULTITHREAD) $(ENABLE_ATOMIC) $(STDF) $(CFLAGS)
HEADER = io.h msa.h

INSTALL = install

MSAMAINOBJ = mtxutl.o io.o constants.o msa_main.o function.o

wmsa : $(MSAMAINOBJ) 
	$(CC) -o $@ $(MSAMAINOBJ) $(MYCFLAGS) $(LIBS)

msa_main.o : msa_main.c $(HEADER)
	$(CC) $(MYCFLAGS) -c msa_main.c 

io.o : io.c $(HEADER)
	$(CC) $(MYCFLAGS) -c io.c 

mtxutl.o : ./submafft/mtxutl.c 
	$(CC) $(MYCFLAGS) -c ./submafft/mtxutl.c 

constants.o : constants.c
	$(CC) $(MYCFLAGS) -c constants.c

function.o : function.c function.h
	$(CC) $(MYCFLAGS) -c function.c

clean:
	rm *.o cd-hit/*.o submafft/*.o tmp* wmsa cd-hit/cd-hit cd-hit/cd-hit-454 cd-hit/cd-hit-est submafft/staralign submafft/disttbfast
	rm -rf swap

all:
	cd submafft && make -j$(THREADS)
	cd cd-hit && make -j$(THREADS) $(NOOPENMPFLAG)
	make -j$(THREADS)

install: 
	mkdir -p $(CDHITFOLDER) $(MAFFTFOLDER)
	cp wmsa $(INSTALLFOLDER)
	cp ./cd-hit/cd-hit ./cd-hit/cd-hit-est $(CDHITFOLDER)
	cp ./submafft/profilealign ./submafft/staralign $(MAFFTFOLDER)

uninstall:
	rm $(INSTALLFOLDER)/wmsa
	rm -rf $(CDHITFOLDER) $(MAFFTFOLDER)