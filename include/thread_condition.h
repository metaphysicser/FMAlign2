#ifndef _CONDITION_H_
#define _CONDITION_H_
// This file is from project https://github.com/malabz/WMSA and has been used with the author's permission.

#if (defined(__linux__) || defined(__APPLE__) || defined(__MINGW32__) || defined(__MINGW64__))
#include <pthread.h>
#else
#endif
// Encapsulates a mutex and condition variable as state
typedef struct condition
{
    pthread_mutex_t pmutex;
    pthread_cond_t pcond;
}condition_t;

// Functions that operate on the state
int condition_init(condition_t* cond);
int condition_lock(condition_t* cond);
int condition_unlock(condition_t* cond);
int condition_wait(condition_t* cond);
int condition_timedwait(condition_t* cond, const struct timespec* abstime);
int condition_signal(condition_t* cond);
int condition_broadcast(condition_t* cond);
int condition_destroy(condition_t* cond);

#endif
