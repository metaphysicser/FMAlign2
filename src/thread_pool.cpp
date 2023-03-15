#include "../include/thread_pool.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <time.h>
#if (defined(__linux__) || defined(__APPLE__))
#include <sys/time.h>
#endif
// This file is from project https://github.com/malabz/WMSA and has been used with the author's permission.

#ifdef _GNU_SOURCE
# undef  _XOPEN_SOURCE
# define _XOPEN_SOURCE 600
# undef  _XOPEN_SOURCE_EXTENDED
# define _XOPEN_SOURCE_EXTENDED 1
# undef  _LARGEFILE64_SOURCE
# define _LARGEFILE64_SOURCE 1
# undef  _BSD_SOURCE
# define _BSD_SOURCE 1
# undef  _SVID_SOURCE
# define _SVID_SOURCE 1
# undef  _ISOC99_SOURCE
# define _ISOC99_SOURCE 1
# undef  _POSIX_SOURCE
# define _POSIX_SOURCE 1
# undef  _POSIX_C_SOURCE
# define _POSIX_C_SOURCE 200112L
# undef  _ATFILE_SOURCE
# define _ATFILE_SOURCE 1
#endif

#if (_WIN32 || _WIN64)
// this block code is written by https://stackoverflow.com/questions/5404277/porting-clock-gettime-to-windows
LARGE_INTEGER getFILETIMEoffset()
{
    SYSTEMTIME s;
    FILETIME f;
    LARGE_INTEGER t;

    s.wYear = 1970;
    s.wMonth = 1;
    s.wDay = 1;
    s.wHour = 0;
    s.wMinute = 0;
    s.wSecond = 0;
    s.wMilliseconds = 0;
    SystemTimeToFileTime(&s, &f);
    t.QuadPart = f.dwHighDateTime;
    t.QuadPart <<= 32;
    t.QuadPart |= f.dwLowDateTime;
    return (t);
}

int clock_gettime(struct timespec* tv)
{
    LARGE_INTEGER           t;
    FILETIME            f;
    double                  microseconds;
    static LARGE_INTEGER    offset;
    static double           frequencyToMicroseconds;
    static int              initialized = 0;
    static BOOL             usePerformanceCounter = 0;

    if (!initialized) {
        LARGE_INTEGER performanceFrequency;
        initialized = 1;
        usePerformanceCounter = QueryPerformanceFrequency(&performanceFrequency);
        if (usePerformanceCounter) {
            QueryPerformanceCounter(&offset);
            frequencyToMicroseconds = (double)performanceFrequency.QuadPart / 1000000.;
        }
        else {
            offset = getFILETIMEoffset();
            frequencyToMicroseconds = 10.;
        }
    }
    if (usePerformanceCounter) QueryPerformanceCounter(&t);
    else {
        GetSystemTimeAsFileTime(&f);
        t.QuadPart = f.dwHighDateTime;
        t.QuadPart <<= 32;
        t.QuadPart |= f.dwLowDateTime;
    }

    t.QuadPart -= offset.QuadPart;
    microseconds = (double)t.QuadPart / frequencyToMicroseconds;
    t.QuadPart = microseconds;
    tv->tv_sec = t.QuadPart / 1000000;
    tv->tv_nsec = t.QuadPart % 1000000;
    return (0);
}
#endif

// The created thread executes
void* thread_routine(void* arg)
{
    struct timespec abstime;
    int timeout;
    // printf("thread %d is starting\n", (int)pthread_self());
    threadpool_t* pool = (threadpool_t*)arg;
    while (1)
    {
        timeout = 0;
        // Need to lock before accessing the thread pool
        condition_lock(&pool->ready);
        // idle
        pool->idle++;
        // Waiting for tasks to arrive in the queue or receiving a thread pool destruction notification
        while (pool->first == NULL && !pool->quit)
        {
            // Otherwise the thread blocks and waits
            // printf("thread %d is waiting\n", (int)pthread_self());
            // Get the current time and add the waiting time to set the timeout sleep time of the process
#if (_WIN32 || _WIN64)
            clock_gettime(&abstime);
#else
            clock_gettime(CLOCK_REALTIME, &abstime);
#endif
            abstime.tv_sec += 10;
            
            int status;
            status = condition_timedwait(&pool->ready, &abstime);  // This function will unlock, allow other threads to access, and when awakened, lock
            if (status == ETIMEDOUT)
            {
//#if (_WIN32 || _WIN64)
//                printf("thread %d wait timed out\n", pthread_self().x);
//#else
#if DEBUG
                printf("thread wait timed out\n");
#endif
// #endif
                timeout = 1;
                break;
            }
        }

        pool->idle--;
        if (pool->first != NULL)
        {
            // Take out the task at the top of the waiting queue, remove the task, and execute the task
            task_t* t = pool->first;
            pool->first = t->next;
            // Since task execution takes time, unlock it first to allow other threads to access the thread pool
            condition_unlock(&pool->ready);
            // perform tasks
            t->run(t->arg);
            // Release the memory after executing the task
            free(t);
            // re-lock
            condition_lock(&pool->ready);
        }

        // Exit the thread pool
        if (pool->quit && pool->first == NULL)
        {
            pool->counter--; // Number of currently working threads - 1
            // If there are no threads in the thread pool, notify the waiting thread (main thread) that all tasks have been completed
            if (pool->counter == 0)
            {
                condition_signal(&pool->ready);
            }
            condition_unlock(&pool->ready);
            break;
        }
        // Timeout, jump out of the destroy thread
        if (timeout == 1)
        {
            pool->counter--;// Number of currently working threads - 1
            condition_unlock(&pool->ready);
            break;
        }

        condition_unlock(&pool->ready);
    }

    // printf("thread %d is exiting\n", (int)pthread_self());
    return NULL;

}


// Thread pool initialization
void threadpool_init(threadpool_t* pool, int threads)
{

    condition_init(&pool->ready);
    pool->first = NULL;
    pool->last = NULL;
    pool->counter = 0;
    pool->idle = 0;
    pool->max_threads = threads;
    pool->quit = 0;

}


// Add a task to the thread pool
void threadpool_add_task(threadpool_t* pool, void* (*run)(void* arg), void* arg)
{
    // spawn a new task
    task_t* newtask = (task_t*)malloc(sizeof(task_t));
    newtask->run = run;
    newtask->arg = arg;
    newtask->next = NULL;// Newly added tasks are placed at the end of the queue

    // The state of the thread pool is shared by multiple threads, and it needs to be locked before the operation
    condition_lock(&pool->ready);

    if (pool->first == NULL)// first task added
    {
        pool->first = newtask;
    }
    else
    {
        pool->last->next = newtask;
    }
    pool->last = newtask;  // The tail of the queue points to the newly joined thread
    // There are idle threads in the thread pool, wake up
    if (pool->idle > 0)
    {
        condition_signal(&pool->ready);
    }
    // The number of threads in the current thread pool does not reach the set maximum value, create a new thread
    else if (pool->counter < pool->max_threads)
    {
        pthread_t tid;
        pthread_create(&tid, NULL, thread_routine, pool);
        pool->counter++;
    }
    
    condition_unlock(&pool->ready);
}

// thread pool destruction
void threadpool_destroy(threadpool_t* pool)
{
    // If destroy has been called, return directly
    if (pool->quit)
    {
        return;
    }
    // lock
    condition_lock(&pool->ready);
    // Set the destroy flag to 1
    pool->quit = 1;
    // The number of threads in the thread pool is greater than 0
    if (pool->counter > 0)
    {
        // For waiting threads, send a signal to wake up
        if (pool->idle > 0)
        {
            condition_broadcast(&pool->ready);
        }
        // Threads that are executing tasks, waiting for them to finish their tasks
        while (pool->counter)
        {
            condition_wait(&pool->ready);
        }
    }
    condition_unlock(&pool->ready);
    condition_destroy(&pool->ready);
}