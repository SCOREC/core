/****************************************************************************** 

  Copyright 2015 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/
#include "reel.h"
#include "reel_config.h"
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <stddef.h>

void reel_fail(const char* format, ...)
{
  fprintf(stderr, "PUMI error: ");
  va_list ap;
  va_start(ap, format);
  vfprintf(stderr, format, ap);
  va_end(ap);
  fprintf(stderr, "\n");
  abort();
}

#if ENABLE_THREADS == 1
#include <pthread.h>

/* Some people don't implement pthread_barrier,
   and some are hard to detect via the C preprocessor,
   so we'll use our own implementation based on other pthread features. */

typedef struct
{
    pthread_mutex_t mutex;
    pthread_cond_t cond;
    int count;
    int tripCount;
} barrier_t;

static void barrier_init(barrier_t *barrier, unsigned int count)
{
  int err;
  err = pthread_mutex_init(&barrier->mutex, 0);
  if (err)
    reel_fail("pthread_mutex_init failed in reel's barrier_init");
  err = pthread_cond_init(&barrier->cond, 0);
  if (err)
    reel_fail("pthread_cond_init failed in reel's barrier_init");
  barrier->tripCount = count;
  barrier->count = 0;
}

static void barrier_destroy(barrier_t *barrier)
{
  pthread_cond_destroy(&barrier->cond);
  pthread_mutex_destroy(&barrier->mutex);
}

static void barrier_wait(barrier_t *barrier)
{
  pthread_mutex_lock(&barrier->mutex);
  ++(barrier->count);
  if (barrier->count > barrier->tripCount)
    reel_fail("bug in reel's barrier_wait");
  if(barrier->count == barrier->tripCount)
  {
    barrier->count = 0;
    pthread_cond_broadcast(&barrier->cond);
    pthread_mutex_unlock(&barrier->mutex);
  }
  else
  {
    pthread_cond_wait(&barrier->cond, &(barrier->mutex));
    pthread_mutex_unlock(&barrier->mutex);
  }
}

#define MAX_THREADS 1024

static pthread_t       global_threads[MAX_THREADS];
static int             global_nthreads = 0;
static pthread_key_t   global_key;
static barrier_t       global_barrier;
static pthread_mutex_t global_lock;

void reel_run_threads(int count, void* (*thread_main)(void*))
{
  if (count < 1)
    reel_fail("reel_run_threads: count must be positive");
  if (count > MAX_THREADS)
    reel_fail("reel_run_threads: count must be less than %d", MAX_THREADS);
  global_nthreads = count;
  global_threads[0] = pthread_self();
  barrier_init(&global_barrier, count);
  pthread_mutex_init(&global_lock, NULL);
  int err;
  err = pthread_key_create(&global_key,NULL);
  if (err)
    reel_fail("pthread_key_create failed");
  pthread_setspecific(global_key,0);
  for (int i = 1; i < count; ++i) {
    err = pthread_create(global_threads + i,
        NULL, thread_main, (void*)(ptrdiff_t)i);
    if (err)
      reel_fail("pthread_create failed");
  }
  thread_main(NULL);
  for (int i = 1; i < count; ++i) {
    err = pthread_join(global_threads[i],NULL);
    if (err)
      reel_fail("pthread_join failed");
  }
  pthread_mutex_destroy(&global_lock);
  barrier_destroy(&global_barrier);
  err = pthread_key_delete(global_key);
  if (err)
    reel_fail("pthread_key_delete failed");
}

int reel_thread_size(void)
{
  if (global_nthreads)
    return global_nthreads;
  return 1;
}

int reel_thread_rank(void)
{
  if (global_nthreads)
    /* double cast to assure GCC that I know what I'm doing. */
    return (int)(ptrdiff_t)(pthread_getspecific(global_key));
  return 0;
}

void reel_thread_init(void* in)
{
  pthread_setspecific(global_key,in);
}

void reel_thread_barrier(void)
{
  if (global_nthreads)
    barrier_wait(&global_barrier);
}

void reel_thread_lock(void)
{
  if (global_nthreads)
    pthread_mutex_lock(&global_lock);
}

void reel_thread_unlock(void)
{
  if (global_nthreads)
    pthread_mutex_unlock(&global_lock);
}
#else
void reel_run_threads(int count, void* (*thread_main)(void*))
{
  (void)count;
  (void)thread_main;
  reel_fail("reel_run_threads: not compiled with -DENABLE_THREADS=ON");
}

int reel_thread_size(void)
{
  return 1;
}

int reel_thread_rank(void)
{
  return 0;
}

void reel_thread_init(void* in)
{
  (void)in;
  reel_fail("reel_thread_init: not compiled with -DENABLE_THREADS=ON");
}

void reel_thread_barrier(void)
{
}

void reel_thread_lock(void)
{
}

void reel_thread_unlock(void)
{
}
#endif


#if defined(__linux__) || defined(__APPLE__)
#include <execinfo.h> /* backtrace for pcu_trace */
#include <signal.h> /* signal for pcu_protect */

void reel_trace(void)
{
  static void* buf[64];
  int n;
  n = backtrace(buf, 64);
  backtrace_symbols_fd(buf, n, 2);
  fflush(stderr);
}

static void catch(int s)
{
  static volatile sig_atomic_t already_dying = 0;
  if (already_dying)
    return;
  already_dying = 1;
  fprintf(stderr, "signal %d caught by pcu\n", s);
  reel_trace();
  signal(s, SIG_DFL);
  raise(s);
}

void reel_protect(void)
{
  signal(SIGABRT, catch);
  signal(SIGSEGV, catch);
  signal(SIGINT, catch);
  signal(SIGFPE, catch);
}
#else
void reel_protect(void)
{
  reel_fail(stderr, "reel_protect only supported on Linux and OS X");
}
#endif
