/****************************************************************************** 

  Copyright 2011 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/
#include <pthread.h>
#include "pcu_thread.h"
#include "pcu_common.h"

static pthread_t* global_threads = NULL;
static int global_nthreads = 0;
static pthread_key_t global_key;
static pthread_barrier_t global_barrier;

void pcu_run_threads(int count, pcu_thread* function)
{
  if (count < 1) pcu_fail("thread count must be positive");
  global_nthreads = count;
  PCU_MALLOC(global_threads,(size_t)count);
  *global_threads = pthread_self();
  pthread_barrier_init(&global_barrier, NULL, count);

  int err;
  err = pthread_key_create(&global_key,NULL);
  if (err) pcu_fail("pthread_key_create failed");
  pthread_setspecific(global_key,0);

  for (int i=1; i < count; ++i)
  {
    err = pthread_create(global_threads+i,NULL,function,NULL);
    if (err) pcu_fail("pthread_create failed");
  }

  function(NULL);
  for (int i=1; i < count; ++i)
  {
    err = pthread_join(global_threads[i],NULL);
    if (err) pcu_fail("pthread_join failed");
  }
  pthread_barrier_destroy(&global_barrier);

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

void pcu_thread_init(void)
{
  int size = pcu_thread_size();
  pthread_t self = pthread_self();
  for (int id=0; id < size; ++id)
  {
    if (pthread_equal(global_threads[id],self))
    {
      pthread_setspecific(global_key,(void*)(ptrdiff_t)id);
      return;
    }
  }
  pcu_fail("could not find thread");
}


void pcu_thread_barrier(void)
{
  pthread_barrier_wait(&global_barrier);
}
