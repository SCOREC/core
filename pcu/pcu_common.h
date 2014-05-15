/****************************************************************************** 

  (c) 2011-2014 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/
#ifndef PCU_COMMON_H
#define PCU_COMMON_H

void pcu_fail(const char* message) __attribute__((noreturn));
int pcu_floor_log2(int n);
int pcu_ceil_log2(int n);

#define MIN(a,b) (((b)<(a))?(b):(a))
#define MAX(a,b) (((b)>(a))?(b):(a))

#endif
