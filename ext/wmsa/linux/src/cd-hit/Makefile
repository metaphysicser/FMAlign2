CC = g++ -Wall -ggdb
CC = g++ -pg
CC = g++

# default with OpenMP
# with OpenMP
# in command line: 
# make openmp=yes
ifeq ($(openmp),no)
  CCFLAGS = -DNO_OPENMP
else
  CCFLAGS = -fopenmp
endif

# prevent `to_string` complie error
CCFLAGS += -std=c++11

#LDFLAGS = -static -lz -o
#LDFLAGS = /usr/lib/x86_64-linux-gnu/libz.a -o

# default with zlib
# without zlib
# in command line:
# make zlib=no
ifeq ($(zlib),no)
  CCFLAGS += 
  LDFLAGS += -o
else
  CCFLAGS += -DWITH_ZLIB
  LDFLAGS += -lz -o
endif

# support debugging
# in command line:
# make debug=yes
# make openmp=yes debug=yes
ifeq ($(debug),yes)
CCFLAGS += -ggdb
else
CCFLAGS += -O2
endif

ifdef MAX_SEQ
CCFLAGS += -DMAX_SEQ=$(MAX_SEQ)
endif

PROGS = cd-hit cd-hit-est cd-hit-454

# Propagate hardening flags
CCFLAGS := $(CPPFLAGS) $(CCFLAGS) $(CXXFLAGS)

.c++.o:
	$(CC) $(CCFLAGS) -c $<

all: $(PROGS)

clean:
	rm -f *.o $(PROGS)

# programs

cd-hit: cdhit-common.o cdhit-utility.o cdhit.o
	$(CC) $(CCFLAGS) cdhit.o cdhit-common.o cdhit-utility.o $(LDFLAGS) cd-hit

cd-hit-est: cdhit-common.o cdhit-utility.o cdhit-est.o
	$(CC) $(CCFLAGS) cdhit-est.o cdhit-common.o cdhit-utility.o $(LDFLAGS) cd-hit-est

cd-hit-454: cdhit-common.o cdhit-utility.o cdhit-454.o
	$(CC) $(CCFLAGS) cdhit-454.o cdhit-common.o cdhit-utility.o $(LDFLAGS) cd-hit-454

# objects
cdhit-common.o: cdhit-common.c++ cdhit-common.h
	$(CC) $(CCFLAGS) cdhit-common.c++ -c

cdhit-utility.o: cdhit-utility.c++ cdhit-utility.h
	$(CC) $(CCFLAGS) cdhit-utility.c++ -c

cdhit.o: cdhit.c++ cdhit-utility.h
	$(CC) $(CCFLAGS) cdhit.c++ -c

cdhit-est.o: cdhit-est.c++ cdhit-utility.h
	$(CC) $(CCFLAGS) cdhit-est.c++ -c

cdhit-454.o: cdhit-454.c++ cdhit-common.h
	$(CC) $(CCFLAGS) cdhit-454.c++ -c


install:
	@echo "only for submodules, not supported install"
