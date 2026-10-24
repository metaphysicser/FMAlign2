cmake_minimum_required(VERSION 3.15)
project(FMAlign2 VERSION 1.0 LANGUAGES CXX C)

# 设置 C++ 标准
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_C_STANDARD 11)

# 编译选项
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall")

# 包含目录
include_directories(include)

# 源文件
set(SOURCES
    main.cpp
    src/utils.cpp
    src/mem_finder.cpp
    src/sequence_split_align.cpp
    src/ssw.cpp
    src/ssw_cpp.cpp
    src/gsacak.c
)

# 根据操作系统添加不同的源文件
if(WIN32)
    # Windows: 使用 OpenMP
    find_package(OpenMP REQUIRED)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
else()
    # Linux/Unix: 使用 pthread 和自定义线程池
    list(APPEND SOURCES
        src/thread_pool.cpp
        src/thread_condition.cpp
    )
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
    if(UNIX AND NOT APPLE)
        # Linux 特定: 需要 rt 库
        set(EXTRA_LIBS rt)
    endif()
endif()

# M64 宏定义 (可选，用于大数据)
option(USE_M64 "Use 64-bit integer types" OFF)
if(USE_M64)
    add_definitions(-DM64=1)
endif()

# Debug 宏定义
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_definitions(-DDEBUG=1)
endif()

# 创建可执行文件
add_executable(fmalign2 ${SOURCES})

# 链接库
if(NOT WIN32)
    target_link_libraries(fmalign2 pthread ${EXTRA_LIBS})
endif()

# 静态链接选项 (Linux)
option(STATIC_LINK "Build with static linking" OFF)
if(STATIC_LINK AND UNIX AND NOT APPLE)
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static -static-libstdc++ -static-libgcc")
    target_link_libraries(fmalign2 -Wl,--whole-archive -lpthread -Wl,--no-whole-archive)
endif()

# Sanitizer 选项 (用于调试)
option(USE_SANITIZER "Build with address and undefined behavior sanitizers" OFF)
if(USE_SANITIZER)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fsanitize=undefined -fno-omit-frame-pointer")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fsanitize=address -fsanitize=undefined -fno-omit-frame-pointer")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fsanitize=address -fsanitize=undefined")
endif()

# 安装配置
install(TARGETS fmalign2 DESTINATION bin)

include(CTest)
if(BUILD_TESTING)
    add_executable(sequence_symbols_test tests/sequence_symbols.cpp src/utils.cpp)
    add_test(NAME sequence_symbols
             COMMAND sequence_symbols_test ${CMAKE_CURRENT_SOURCE_DIR}/tests/symbols.fasta)
    set(PIPELINE_SOURCES ${SOURCES})
    list(REMOVE_ITEM PIPELINE_SOURCES main.cpp)
    add_executable(symbols_pipeline_test tests/symbols_pipeline.cpp ${PIPELINE_SOURCES})
    if(NOT WIN32)
        target_link_libraries(symbols_pipeline_test pthread ${EXTRA_LIBS})
    endif()

    # This test includes mem_finder.cpp with a test-only allocation shim.
    add_executable(lcp_initialization_test
        tests/lcp_initialization.cpp src/utils.cpp src/gsacak.c)
    if(NOT WIN32)
        target_sources(lcp_initialization_test PRIVATE
            src/thread_pool.cpp src/thread_condition.cpp)
        target_link_libraries(lcp_initialization_test pthread ${EXTRA_LIBS})
    endif()
    add_test(NAME lcp_initialization COMMAND lcp_initialization_test)

    find_package(Python3 COMPONENTS Interpreter QUIET)
    if(Python3_Interpreter_FOUND)
        add_test(NAME symbols_pipeline
            COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_symbols_cli.py
                    $<TARGET_FILE:symbols_pipeline_test> --pipeline)
        add_test(NAME symbols_cli
            COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_symbols_cli.py
                    $<TARGET_FILE:fmalign2>)
        set_tests_properties(symbols_pipeline symbols_cli PROPERTIES TIMEOUT 180)
    endif()
endif()

