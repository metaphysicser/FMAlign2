CC = g++
CFLAGS = -Wall -g

SRCS = main.cpp src\utils.cpp src\mem_finder.cpp src\mem_filter.cpp src\sequence_split_align.cpp
OBJS = $(SRCS:.cpp=.o)

FMAlign2: $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o FMAlign2

main.o: main.cpp include\utils.h include\common.h include\mem_finder.h include\mem_filter.h include\sequence_split_align.h
	$(CC) $(CFLAGS) -c main.cpp

utils.o: utils.cpp include\utils.h include\common.h include\kseq.h
	$(CC) $(CFLAGS) -c utils.cpp

mem_finder.o: mem_finder.cpp include\common.h
	$(CC) $(CFLAGS) -c mem_finder.cpp

mem_filter.o: mem_filter.cpp include\common.h
	$(CC) $(CFLAGS) -c mem_filter.cpp

sequence_split_align.o: sequence_split_align.cpp include\common.h
	$(CC) $(CFLAGS) -c sequence_split_align.cpp
clean:
	
ifeq ($(OS),Windows_NT)
	del -f $(OBJS) FMAlign2.exe
else
	rm -f $(OBJS) FMAlign2
endif
	
