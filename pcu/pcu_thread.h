/****************************************************************************** 

  (c) 2011-2014 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/
#ifndef PCU_THREAD_H
#define PCU_THREAD_H

#include "pcu_memory.h"

typedef void* pcu_thread(void*);

void pcu_run_threads(int count, pcu_thread* function);
int pcu_thread_size(void);
int pcu_thread_rank(void);
void pcu_thread_init(void);
void pcu_thread_barrier (void);

#endif
