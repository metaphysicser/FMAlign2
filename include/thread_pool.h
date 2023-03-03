#ifndef _THREAD_POOL_H_
#define _THREAD_POOL_H_
// This file is from project https://github.com/malabz/WMSA and has been used with the author's permission.
// Thread pool header file
#if (_WIN32 || _WIN64)
#include <Windows.h>
#endif

#include "thread_condition.h"

// Encapsulate the task object that needs to be executed by the object in the thread pool
typedef struct task
{
    void* (*run)(void* args);  // Function pointer, the task that needs to be performed
    void* arg;              // parameter
    struct task* next;      // The next task in the task queue
}task_t;


//下面是线程池结构体
typedef struct threadpool
{
    condition_t ready;    // state quantity
    task_t* first;       // The first task in the task queue
    task_t* last;        // The last task in the task queue
    int counter;         // The number of existing threads in the thread pool
    int idle;            // The number of idle threads in the thread pool
    int max_threads;     // The maximum number of threads in the thread pool
    int quit;            //Whether to exit the flag
}threadpool_t;


// Thread pool initialization
void threadpool_init(threadpool_t* pool, int threads);

// Add tasks to the thread pool
void threadpool_add_task(threadpool_t* pool, void* (*run)(void* arg), void* arg);

// destroy thread pool
void threadpool_destroy(threadpool_t* pool);

#endif
