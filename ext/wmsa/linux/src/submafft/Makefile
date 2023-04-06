PREFIX = /usr/local
LIBDIR = $(PREFIX)/libexec/mafft
BINDIR = $(PREFIX)/bin
MANDIR = $(PREFIX)/share/man/man1
DESTDIR = 

#MNO_CYGWIN = -mno-cygwin

# ENABLE_MULTITHREAD = -Denablemultithread
# Comment out the above line if your compiler 
# does not support TLS (thread-local strage).

#ENABLE_ATOMIC = -Denableatomic
# Comment out the above line if your compiler 
# does not support "atomic_int".

#DASH_CLIENT = dash_client
# Uncomment the above line to use protein 3D 
# structural information.  Go language is required.

CC = gcc
#CC = icc
CFLAGS = -g -O3 -D _GNU_SOURCE -m64
#CFLAGS = -O3 -fPIC
# add -fPIC when building .so files

#CC = icc
#CFLAGS = -fast
# if you have icc, use this.

#CFLAGS =  -O0  -fPIC -pedantic -Wall -std=c99 -g -pg -DMALLOC_CHECK_=3
#CFLAGS =  -O0  -fPIC -pedantic -Wall -std=c99 -g -pg -DMALLOC_CHECK_=3  -fprofile-arcs -ftest-coverage 
#CFLAGS =  -O0  -fPIC -pedantic -Wall -std=c99 -g -DMALLOC_CHECK_=3 # for shark, valgrind
#CFLAGS =  -O0  -fPIC -pedantic -Wall -std=c99 -g -DMALLOC_CHECK_=3 -lprofiler  # ? 

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

INSTALL = install

STRIP = strip
#STRIP = true # to disable strip

PROGS = staralign profilealign 
SOS = libdisttbfast.so
DLLS = libdisttbfast.dll
DYLIBS = libdisttbfast.dylib

SCRIPTS = 
OBJSTARALIGN = mtxutl.o io.o mltaln9.o tddis.o constants.o \
		    Falign.o Galign11.o SAalignmm.o \
			staralign.o defs.o fft.o fftFunctions.o  \
			 Salignmm.o Kband.o $(MULTIOBJ) 
OBJPROFILEALIGN = mtxutl.o io.o mltaln9.o tddis.o constants.o \
		    Falign.o Galign11.o SAalignmm.o \
			profilealign.o defs.o fft.o fftFunctions.o  \
			 Salignmm.o Kband.o $(MULTIOBJ)
ifdef ENABLE_MULTITHREAD
MULTIOBJ = threadpool.o threadpool_condition.o
endif

HEADER = mltaln.h mtxutl.h 
ifdef ENABLE_MULTITHREAD
MULTIHEADER = threadpool.h threadpool_condition.h
endif
FFTHEADER = fft.h
KBANDHEADER = Kband.h


all : $(PROGS) $(SCRIPTS)
	@echo done.

sos : $(SOS)
dylibs : $(DYLIBS)
dlls : $(DLLS)

$(DASH_CLIENT): dash_client.go
	go build dash_client.go
#	go build --ldflags '-extldflags "-static"' dash_client.go

univscript: univscript.tmpl Makefile
	sed "s:_PROGS:$(PROGS):" univscript.tmpl  > univscript

mafft: mafft.tmpl mltaln.h
	sed "s:_LIBDIR:$(LIBDIR):" mafft.tmpl  > mafft

mltaln.h : functions.h
	touch mltaln.h

version : version.c mltaln.h
	$(CC) -o $@ version.c $(MYCFLAGS) $(LDFLAGS) $(LIBS) 


staralign: $(OBJSTARALIGN)
	$(CC) -o $@ $(OBJSTARALIGN) $(MYCFLAGS) $(LDFLAGS) $(LIBS)

profilealign: $(OBJPROFILEALIGN)
	$(CC) -o $@ $(OBJPROFILEALIGN) $(MYCFLAGS) $(LDFLAGS) $(LIBS)

libdisttbfast.so : $(OBJDISTTBFAST)
	$(CC) -shared -o $@ $(OBJDISTTBFAST) $(MYCFLAGS) $(LDFLAGS) $(LIBS)

libdisttbfast.dylib : $(OBJDISTTBFAST)
	$(CC) -dynamiclib -o $@ $(OBJDISTTBFAST) $(MYCFLAGS) $(LDFLAGS) $(LIBS)

libdisttbfast.dll : $(OBJDISTTBFAST)
	$(CC) -shared -o $@ $(OBJDISTTBFAST) $(MYCFLAGS) $(LDFLAGS) $(LIBS)


mltaln9.o : mltaln9.c $(HEADER)
	$(CC) $(MYCFLAGS) -c mltaln9.c

tddis.o : tddis.c $(HEADER)
	$(CC) $(MYCFLAGS) -c tddis.c

constants.o : constants.c miyata.h miyata5.h blosum.c DNA.h JTT.c $(HEADER)
	$(CC) $(MYCFLAGS) -c constants.c

defs.o : defs.c 
	$(CC) $(MYCFLAGS) -c defs.c

Salignmm.o : Salignmm.c $(HEADER)
	$(CC) $(MYCFLAGS) -c Salignmm.c 

Galign11.o : Galign11.c $(HEADER)
	$(CC) $(MYCFLAGS) -c Galign11.c 

SAalignmm.o : SAalignmm.c $(HEADER)
	$(CC) $(MYCFLAGS) -c SAalignmm.c -o SAalignmm.o

disttbfast.o : disttbfast.c $(HEADER) $(FFTHEADER)
	$(CC) $(MYCFLAGS) -c disttbfast.c

staralign.o : staralign.c $(HEADER) $(FFTHEADER) $(MULTIHEADER)
	$(CC) $(MYCFLAGS) -c staralign.c

profilealign.o : profilealign.c $(HEADER) $(FFTHEADER) $(MULTIHEADER)
	$(CC) $(MYCFLAGS) -c profilealign.c

threadpool_condition.o : threadpool_condition.c $(HEADER) $(MULTIHEADER)
	$(CC) $(MYCFLAGS) -c threadpool_condition.c

threadpool.o : threadpool.c $(HEADER) $(MULTIHEADER)
	$(CC) $(MYCFLAGS) -c threadpool.c

io.o : io.c $(HEADER) $(FFTHEADER)
	$(CC) $(MYCFLAGS) -c io.c

fft.o : fft.c $(HEADER) $(FFTHEADER)
	$(CC) $(MYCFLAGS) -c fft.c 

fftFunctions.o : fftFunctions.c $(HEADER) $(FFTHEADER)
	$(CC) $(MYCFLAGS) -c fftFunctions.c

Falign.o : Falign.c $(HEADER) $(FFTHEADER) $(MTXHEADER)
	$(CC) $(MYCFLAGS) -c Falign.c

mtxutl.o : mtxutl.c 
	$(CC) $(MYCFLAGS) -c mtxutl.c

addfunctions.o : addfunctions.c $(HEADER)
	$(CC) $(MYCFLAGS) -c addfunctions.c

Kband.o : Kband.c $(KBANDHEADER)
	$(CC) $(MYCFLAGS) -c Kband.c

clean :
	rm -f *.o *.a *.exe *~ $(PROGS) $(SOS) $(DYLIBS) $(DLLS) *.gcda *.gcno $(DASH_CLIENT) pre trace
#	rm -f ../binaries/* ../scripts/*

install : all
	mkdir -p $(DESTDIR)$(LIBDIR)
	chmod 755 $(DESTDIR)$(LIBDIR)
	mkdir -p $(DESTDIR)$(BINDIR)
	chmod 755 $(DESTDIR)$(BINDIR)
	chmod 755 $(SCRIPTS)
	$(INSTALL) $(SCRIPTS)  $(DESTDIR)$(BINDIR)
	chmod 755 $(PROGS) ||:     # in MinGW, it's ok if this fails
	$(STRIP) $(PROGS) ||: # may fail for dash_client on mac.
	$(INSTALL) $(PROGS) $(DESTDIR)$(LIBDIR)

	( cd $(DESTDIR)$(BINDIR); \
rm -f fftns fftnsi; \
rm -f mafft-fftns mafft-fftnsi; \
ln -s mafft fftns; ln -s mafft fftnsi; \
ln -s mafft mafft-fftns; \
ln -s mafft mafft-fftnsi;)

	mkdir -p $(DESTDIR)$(MANDIR)
	chmod 755 $(DESTDIR)$(MANDIR) 
