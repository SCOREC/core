/****************************************************************************** 

  Copyright 2011 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/
#include <pthread.h>
#include "pcu_thread.h"
#include "pcu_common.h"

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

static void barrier_init(barrier_t *barrier, unsigned int count);
static void barrier_destroy(barrier_t *barrier);
static void barrier_wait(barrier_t *barrier);

static pthread_t* global_threads = NULL;
static int global_nthreads = 0;
static pthread_key_t global_key;
static barrier_t global_barrier;
static pthread_mutex_t global_lock;

void pcu_run_threads(int count, pcu_thread* function)
{
  if (count < 1) pcu_fail("thread count must be positive");
  global_nthreads = count;
  PCU_MALLOC(global_threads,(size_t)count);
  *global_threads = pthread_self();
  barrier_init(&global_barrier, count);
  pthread_mutex_init(&global_lock, NULL);

  int err;
  err = pthread_key_create(&global_key,NULL);
  if (err) pcu_fail("pthread_key_create failed");
  pthread_setspecific(global_key,0);

  for (int i=1; i < count; ++i)
  {
    err = pthread_create(global_threads+i,NULL,function,(void*)(ptrdiff_t)i);
    if (err) pcu_fail("pthread_create failed");
  }

  function(NULL);
  for (int i=1; i < count; ++i)
  {
    err = pthread_join(global_threads[i],NULL);
    if (err) pcu_fail("pthread_join failed");
  }
  pthread_mutex_destroy(&global_lock);
  barrier_destroy(&global_barrier);

  err = pthread_key_delete(global_key);
  if (err) pcu_fail("pthread_key_delete failed");
  pcu_free(global_threads);
}

int pcu_thread_size(void)
{
  return global_nthreads;
}

int pcu_thread_rank(void)
{
  /* double cast to assure GCC that I know what I'm doing. */
  return (int)(ptrdiff_t)(pthread_getspecific(global_key));
}

void pcu_thread_init(void* in)
{
  pthread_setspecific(global_key,in);
}

void pcu_thread_barrier(void)
{
  barrier_wait(&global_barrier);
}

void pcu_thread_lock(void)
{
  pthread_mutex_lock(&global_lock);
}

void pcu_thread_unlock(void)
{
  pthread_mutex_unlock(&global_lock);
}

static void barrier_init(barrier_t *barrier, unsigned int count)
{
  int err;
  err = pthread_mutex_init(&barrier->mutex, 0);
  if (err)
    pcu_fail("mutex_init failed in barrier_init");
  err = pthread_cond_init(&barrier->cond, 0);
  if (err)
    pcu_fail("cond_init failed in barrier_init");
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
    pcu_fail("BUG");
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
