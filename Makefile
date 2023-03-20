CXX = g++
CXXFLAGS = -Wall -std=c++11
CC = g++
CFLAGS = 

SRCS = main.cpp src/utils.cpp src/mem_finder.cpp src/sequence_split_align.cpp ext/SW/ssw.cpp ext/SW/ssw_cpp.cpp

ifdef DEBUG
	CXXFLAGS += -O0 -g -DDEBUG
	CFLAGS += -O0 -g
else
	CXXFLAGS += -O2
	CFLAGS += -O2
endif

ifeq ($(OS),Windows_NT)
    CXXFLAGS +=  -fopenmp
else
	CXXFLAGS += -lpthread -pthread -lrt
	SRCS += src/thread_pool.cpp src/thread_condition.cpp
endif
OBJS = $(SRCS:.cpp=.o)
OBJS += src/gsacak.o

ifdef M64
	CXXFLAGS += -DM64
	CFLAGS += -DM64
endif

FMAlign2: $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o FMAlign2

utils.o: src/utils.cpp include/utils.h include/common.h include/kseq.h
	$(CXX) $(CXXFLAGS) -c src/utils.cpp -o $@

mem_finder.o: src/mem_finder.cpp include/common.h include/gsacak.h include/thread_pool.h
	$(CXX) $(CXXFLAGS) -c src/mem_finder.cpp -o $@

sequence_split_align.o: src/sequence_split_align.cpp include/common.h
	$(CXX) $(CXXFLAGS) -c src/sequence_split_align.cpp -o $@

gsacak.o: src/gsacak.c include/gsacak.h
	$(C) $(CFLAGS)  -c src/gsacak.c -o $@

ifeq ($(OS),Windows_NT)
else
thread_pool.o: src/thread_pool.cpp include/thread_pool.h include/thread_condition.h
	$(CXX) $(CXXFLAGS) -c src/thread_pool.cpp -o $@

thread_condition.o: src/thread_condition.cpp include/thread_condition.h
	$(CXX) $(CXXFLAGS) -c src/thread_conition.cpp -o $@
endif

ssw.o: ext/SW/ssw.cpp ext/SW/ssw.h
	$(CXX) $(CXXFLAGS) -c ext/SW/ssw.cpp -o $@

ssw_cpp.o: ext/SW/ssw_cpp.cpp ext/SW/ssw_cpp.h ext/SW/ssw.h
	$(CXX) $(CXXFLAGS) -c ext/SW/ssw_cpp.cpp -o $@
	
main.o: main.cpp include/utils.h include/common.h include/mem_finder.h include/sequence_split_align.h
	$(CXX) $(CXXFLAGS) -c main.cpp -o $@

clean:
ifeq ($(OS),Windows_NT)
	del -f $(subst /,\\,$(OBJS)) FMAlign2.exe
else
	rm -f $(OBJS) FMAlign2
endif
	
