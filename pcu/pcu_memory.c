/****************************************************************************** 

  Copyright 2011 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/
#include <stdlib.h>
#include "pcu_common.h"
#include "pcu_memory.h"
#include "pcu_thread.h"

//#define DEBUG_MEMORY

#ifdef DEBUG_MEMORY
#include <pthread.h>
static size_t allocated = 0;
pthread_mutex_t allocated_mutex = PTHREAD_MUTEX_INITIALIZER;
#endif

void* pcu_malloc(size_t size)
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
  if (!p) pcu_fail("malloc ran out of memory");
#ifdef DEBUG_MEMORY
  *p = size;
  ++p;
#endif
  return p;
}

void* pcu_realloc(void* p, size_t size)
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
  p2 = realloc(p2,size);
  if ((!p2)&&(size)) pcu_fail("realloc ran out of memory");
#ifdef DEBUG_MEMORY
  if (p2)
  {
    *p2 = size;
    ++p2;
  }
#endif
  return p2;
}

void pcu_free(void* p)
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

size_t pcu_memory(void)
{
#ifdef DEBUG_MEMORY
  return allocated;
#else
  return 0;
#endif
}

void pcu_make_buffer(pcu_buffer* b)
{
  b->start = NULL;
  b->size = 0;
  b->capacity = 0;
}

void pcu_free_buffer(pcu_buffer* b)
{
  pcu_free(b->start);
}

void* pcu_push_buffer(pcu_buffer* b, size_t size)
{
  b->size += size;
  if (b->size > b->capacity)
  {
    //this growth formula is from git's cache.h alloc_nr
    size_t min_growth = ((b->capacity + 16)*3)/2;
    b->capacity = MAX(b->size,min_growth);
    b->start = pcu_realloc(b->start,b->capacity);
  }
  return b->start + b->size - size;
}

void pcu_begin_buffer(pcu_buffer* b)
{
  b->size = 0;
}

void* pcu_walk_buffer(pcu_buffer* b, size_t size)
{
  void* at = b->start + b->size;
  b->size += size;
  if (b->size > b->capacity) pcu_fail("walked past end of buffer");
  return at;
}

bool pcu_buffer_walked(pcu_buffer* b)
{
  return b->size == b->capacity;
}

void pcu_resize_buffer(pcu_buffer* b, size_t size)
{
  if (b->size == size && b->capacity == size) return;
  b->size = b->capacity = size;
  b->start = pcu_realloc(b->start,size);
}

void pcu_set_buffer(pcu_buffer* b, void* p, size_t size)
{
  b->start = p;
  b->size = b->capacity = size;
}

