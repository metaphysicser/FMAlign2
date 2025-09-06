# compilers
CXX     = g++
CC      = gcc

# flags
CXXFLAGS = -Wall -std=c++17 -Iinclude -static
CFLAGS   = -Iinclude
LDFLAGS  =
LDLIBS   =

SRCS = main.cpp \
       src/utils.cpp \
       src/mem_finder.cpp \
       src/sequence_split_align.cpp \
       src/ssw.cpp \
       src/ssw_cpp.cpp

ifeq ($(STATIC_LINK), 1)
	LDFLAGS = -Wl,--whole-archive -lpthread -Wl,--no-whole-archive
else
	LDFLAGS = -pthread
endif

# non-Windows extra sources
ifeq ($(OS),Windows_NT)
    # OpenMP（按需）
    CXXFLAGS += -fopenmp
else
    SRCS += src/thread_pool.cpp src/thread_condition.cpp
    # POSIX 线程
    CXXFLAGS += -pthread
    LDFLAGS  += -pthread
    # Linux 专用的 rt 库（macOS 没有）
    UNAME_S := $(shell uname -s)
    ifeq ($(UNAME_S),Linux)
        LDLIBS += -lrt
    endif
endif

# C 源单独对象（不要混进 SRCS 以免用错编译器）
CSRCS = src/gsacak.c

OBJS  = $(SRCS:.cpp=.o) $(CSRCS:.c=.o)

# debug / release
ifdef DEBUG
    CXXFLAGS += -O0 -g -DDEBUG
    CFLAGS   += -O0 -g
else
    CXXFLAGS += -O3
    CFLAGS   += -O3
endif

# 64-bit 宏
ifdef M64
    CXXFLAGS += -DM64
    CFLAGS   += -DM64
endif

# final target
fmalign2: $(OBJS)
	$(CXX) $(LDFLAGS) $(OBJS) $(LDLIBS) -o $@

# generic rules（自动处理带路径的 .o）
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# 依赖（可选：生成 .d 依赖）
# CXXFLAGS += -MMD -MP
# CFLAGS   += -MMD -MP
# -include $(OBJS:.o=.d)

# clean
.PHONY: clean
clean:
ifeq ($(OS),Windows_NT)
	- del /f $(subst /,\\,$(OBJS)) 2> NUL
	- del /f fmalign2.exe 2> NUL
else
	rm -f $(OBJS) fmalign2
endif
