#include "../include/thread_condition.h"
// This file is from project https://github.com/malabz/WMSA and has been used with the author's permission.
// initialization
int condition_init(condition_t* cond)
{
    int status;
    if ((status = pthread_mutex_init(&cond->pmutex, NULL)))
        return status;

    if ((status = pthread_cond_init(&cond->pcond, NULL)))
        return status;

    return 0;
}

// lock
int condition_lock(condition_t* cond)
{
    return pthread_mutex_lock(&cond->pmutex);
}

// unlock
int condition_unlock(condition_t* cond)
{
    return pthread_mutex_unlock(&cond->pmutex);
}

// wait
int condition_wait(condition_t* cond)
{
    return pthread_cond_wait(&cond->pcond, &cond->pmutex);
}

// fixed time wait
int condition_timedwait(condition_t* cond, const struct timespec* abstime)
{
    return pthread_cond_timedwait(&cond->pcond, &cond->pmutex, abstime);
}

// wake up a sleeping thread
int condition_signal(condition_t* cond)
{
    return pthread_cond_signal(&cond->pcond);
}

// wake up all sleeping threads
int condition_broadcast(condition_t* cond)
{
    return pthread_cond_broadcast(&cond->pcond);
}

// free
int condition_destroy(condition_t* cond)
{
    int status;
    if ((status = pthread_mutex_destroy(&cond->pmutex)))
        return status;

    if ((status = pthread_cond_destroy(&cond->pcond)))
        return status;

    return 0;
}