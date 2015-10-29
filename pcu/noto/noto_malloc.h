/****************************************************************************** 

  Copyright 2015 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/
#ifndef NOTO_MALLOC_H
#define NOTO_MALLOC_H

#include <stddef.h>

void* noto_malloc(size_t size);
#define NOTO_MALLOC(p,c) ((p)=noto_malloc(sizeof(*(p))*(c)))
void* noto_realloc(void* p, size_t size);
void noto_free(void* p);
size_t noto_malloced(void);

#endif

