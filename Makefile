# compilers
CXX     = g++
CC      = gcc

# toggle: 1 = static (default), 0 = dynamic
STATIC_LINK ?= 1

# detect OS
UNAME_S := $(shell uname -s)

# flags (编译期)
CXXFLAGS = -Wall -std=c++17 -Iinclude
CFLAGS   = -Iinclude

# flags (链接期)
LDFLAGS  =
LDLIBS   =

SRCS = main.cpp \
       src/utils.cpp \
       src/mem_finder.cpp \
       src/sequence_split_align.cpp \
       src/ssw.cpp \
       src/ssw_cpp.cpp

# non-Windows extra sources
ifeq ($(OS),Windows_NT)
    # OpenMP（按需）
    CXXFLAGS += -fopenmp
else
    SRCS += src/thread_pool.cpp src/thread_condition.cpp
    # POSIX 线程
    CXXFLAGS += -pthread
    # Linux 专用的 rt 库（macOS 没有）
    ifeq ($(UNAME_S),Linux)
        LDLIBS += -lrt
    endif
endif

# 链接模式：默认静态
ifeq ($(STATIC_LINK),1)
  ifeq ($(UNAME_S),Linux)
    # 静态链接（glibc 环境需要对应的静态库已安装）
    LDFLAGS += -static -static-libstdc++ -static-libgcc
    # pthread 静态：确保被完整吸入
    LDLIBS  += -Wl,--whole-archive -lpthread -Wl,--no-whole-archive
  else
    $(warning Static linking is not supported on $(UNAME_S); falling back to dynamic)
    LDFLAGS += -pthread
  endif
else
  # 动态链接
  LDFLAGS += -pthread
endif

# C 源单独对象
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

# default target
all: fmalign2

# final target
fmalign2: $(OBJS)
	$(CXX) $(LDFLAGS) $(OBJS) $(LDLIBS) -o $@

# generic rules
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# clean
.PHONY: clean
clean:
ifeq ($(OS),Windows_NT)
	- del /f $(subst /,\\,$(OBJS)) 2> NUL
	- del /f fmalign2.exe 2> NUL
else
	rm -f $(OBJS) fmalign2
endif
