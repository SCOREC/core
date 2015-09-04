/****************************************************************************** 

  Copyright 2011 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/
#include <stdlib.h>
#include "pcu_buffer.h"
#include "noto_malloc.h"
#include "reel.h"

void pcu_make_buffer(pcu_buffer* b)
{
  b->start = NULL;
  b->size = 0;
  b->capacity = 0;
}

void pcu_free_buffer(pcu_buffer* b)
{
  noto_free(b->start);
}

void* pcu_push_buffer(pcu_buffer* b, size_t size)
{
  b->size += size;
  if (b->size > b->capacity)
  {
    //this growth formula is from git's cache.h alloc_nr
    size_t min_growth = ((b->capacity + 16)*3)/2;
    b->capacity = b->size;
    if (min_growth > b->capacity)
      b->capacity = min_growth;
    b->start = noto_realloc(b->start, b->capacity);
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
  if (b->size > b->capacity)
    reel_fail("pcu_walk_buffer: walked past end of buffer");
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
  b->start = noto_realloc(b->start,size);
}

void pcu_set_buffer(pcu_buffer* b, void* p, size_t size)
{
  b->start = p;
  b->size = b->capacity = size;
}

