CXX = g++
CXXFLAGS = -Wall -std=c++11

ifdef DEBUG
	CXXFLAGS += -O0 -g -DDEBUG
endif

ifdef M64
	CXXFLAGS += -DM64
endif

SRCS = main.cpp src\utils.cpp src\mem_finder.cpp src\mem_filter.cpp src\sequence_split_align.cpp
OBJS = $(SRCS:.cpp=.o)

FMAlign2: $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o FMAlign2


utils.o: src\utils.cpp include\utils.h include\common.h include\kseq.h
	$(CXX) $(CXXFLAGS) -c src\utils.cpp -g

mem_finder.o: src\mem_finder.cpp include\common.h
	$(CXX) $(CXXFLAGS) -c src\mem_finder.cpp -g

mem_filter.o: src\mem_filter.cpp include\common.h
	$(CXX) $(CXXFLAGS) -c src\mem_filter.cpp -g

sequence_split_align.o: src\sequence_split_align.cpp include\common.h
	$(CXX) $(CXXFLAGS) -c src\sequence_split_align.cpp -g
	
main.o: main.cpp include\utils.h include\common.h include\mem_finder.h include\mem_filter.h include\sequence_split_align.h
	$(CXX) $(CXXFLAGS) -c main.cpp

clean:
ifeq ($(OS),Windows_NT)
	del -f $(OBJS) FMAlign2.exe
else
	rm -f $(OBJS) FMAlign2
endif
	
