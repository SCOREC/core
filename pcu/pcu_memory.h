/****************************************************************************** 

  (c) 2011-2014 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/
#ifndef PCU_MEMORY_H
#define PCU_MEMORY_H

#include <stddef.h>
#include <stdbool.h>

void* pcu_malloc(size_t size);
#define PCU_MALLOC(p,c) ((p)=pcu_malloc(sizeof(*(p))*(c)))
void* pcu_realloc(void* p, size_t size);
void pcu_free(void* p);
size_t pcu_memory(void);

typedef struct
{
  char* start;
  size_t size;
  size_t capacity;
} pcu_buffer;

void pcu_make_buffer(pcu_buffer* b);
void pcu_free_buffer(pcu_buffer* b);
void* pcu_push_buffer(pcu_buffer* b, size_t size);
void pcu_begin_buffer(pcu_buffer* b);
void* pcu_walk_buffer(pcu_buffer* b, size_t size);
bool pcu_buffer_walked(pcu_buffer* b);
void pcu_resize_buffer(pcu_buffer* b, size_t size);
void pcu_set_buffer(pcu_buffer* b, void* p, size_t size);

#endif
