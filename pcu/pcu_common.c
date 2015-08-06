/****************************************************************************** 

  Copyright 2011 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include "pcu_common.h"

void pcu_fail(const char* format, ...)
{
  fprintf(stderr, "PCU error: ");
  va_list ap;
  va_start(ap, format);
  vfprintf(stderr, format, ap);
  va_end(ap);
  fprintf(stderr, "\n");
  abort();
}

int pcu_floor_log2(int n)
{
  int r = 0;
  while ((n >>= 1)) ++r;
  return r;
}

int pcu_ceil_log2(int n)
{
  int r = pcu_floor_log2(n);
  if ((1 << r)<n) ++r;
  return r;
}

