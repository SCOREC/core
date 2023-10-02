/****************************************************************************** 

  Copyright 2015 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/
#include "noto_malloc.h"
#include "reel.h" /* for reel_fail */
#include <stdlib.h>

//#define DEBUG_MEMORY

#ifdef DEBUG_MEMORY
#include <pthread.h>
static size_t allocated = 0;
static pthread_mutex_t allocated_mutex = PTHREAD_MUTEX_INITIALIZER;
#endif

void* noto_malloc(size_t size)
{
  if (!size)
    return NULL;
#ifdef DEBUG_MEMORY
  size += sizeof(size_t);
  pthread_mutex_lock(&allocated_mutex);
  allocated += size;
  pthread_mutex_unlock(&allocated_mutex);
#endif
  size_t* p = malloc(size);
  if (!p)
    reel_fail("malloc(%lu) failed", (unsigned long) size);
#ifdef DEBUG_MEMORY
  *p = size;
  ++p;
#endif
  return p;
}

void* noto_realloc(void* p, size_t size)
{
  if ((!p)&&(!size))
    return NULL;
  size_t* p2 = p;
#ifdef DEBUG_MEMORY
  size += sizeof(size_t);
  pthread_mutex_lock(&allocated_mutex);
  if (p2)
  {
    --p2;
    allocated -= *p2;
  }
  allocated += size;
  pthread_mutex_unlock(&allocated_mutex);
#endif
  size_t* ret = realloc(p2, size);
  if ((!ret) && (size))
    reel_fail("realloc(%p, %lu) failed", (void*) p2, size);
  p2 = ret;
#ifdef DEBUG_MEMORY
  if (p2)
  {
    *p2 = size;
    ++p2;
  }
#endif
  return p2;
}

void noto_free(void* p)
{
  size_t* p2 = p;
#ifdef DEBUG_MEMORY
  if (p)
  {
    --p2;
    pthread_mutex_lock(&allocated_mutex);
    allocated -= *p2;
    pthread_mutex_unlock(&allocated_mutex);
  }
#endif
  free(p2);
}

size_t noto_malloced(void)
{
#ifdef DEBUG_MEMORY
  pthread_mutex_lock(&allocated_mutex);
  size_t sample = allocated;
  pthread_mutex_unlock(&allocated_mutex);
  return sample;
#else
  return 0;
#endif
}
