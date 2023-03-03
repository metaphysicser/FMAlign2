CXX = g++
CXXFLAGS = -Wall -std=c++11
CC = g++
CFLAGS = 
ifdef DEBUG
	CXXFLAGS += -O0 -g -DDEBUG
	CFLAGS += -O0 -g -DDEBUG
else
	CXXFLAGS += -O3
	CFLAGS += -O3
endif


ifeq ($(OS),Windows_NT)
    CXXFLAGS += -fopenmp -L.\ext\pthread-win32 -lpthread
else
  
	CXXFLAGS += -fopenmp

endif


SRCS = main.cpp src/utils.cpp src/mem_finder.cpp src/mem_filter.cpp src/sequence_split_align.cpp src/thread_pool.cpp src/thread_condition.cpp
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

mem_finder.o: src/mem_finder.cpp include/common.h include/gsacak.h
	$(CXX) $(CXXFLAGS) -c src/mem_finder.cpp -o $@

mem_filter.o: src/mem_filter.cpp include/common.h
	$(CXX) $(CXXFLAGS) -c src/mem_filter.cpp -o $@

sequence_split_align.o: src/sequence_split_align.cpp include/common.h
	$(CXX) $(CXXFLAGS) -c src/sequence_split_align.cpp -o $@

gsacak.o: src/gsacak.c include/gsacak.h
	$(C) $(CFLAGS)  -c src/gsacak.c -o $@

thread_pool.o: src/thread_pool.cpp include/thread_pool.h include/thread_condition.h
	$(CXX) $(CXXFLAGS) -c src/thread_pool.cpp -o $@

thread_condition.o: src/thread_condition.cpp include/thread_condition.h
	$(CXX) $(CXXFLAGS) -c src/thread_conition.cpp -o $@
	
main.o: main.cpp include/utils.h include/common.h include/mem_finder.h include/mem_filter.h include/sequence_split_align.h
	$(CXX) $(CXXFLAGS) -c main.cpp -o $@

clean:
ifeq ($(OS),Windows_NT)
	del -f $(subst /,\\,$(OBJS)) FMAlign2.exe
else
	rm -f $(SRCS:.cpp=.o) FMAlign2
endif
	
