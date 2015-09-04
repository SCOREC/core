/****************************************************************************** 

  Copyright 2015 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/
#ifndef REEL_H
#define REEL_H

void reel_fail(const char* format, ...)
  __attribute__((noreturn,format(printf,1,2)));

void reel_run_threads(int count, void* (*thread_main)(void*));
int reel_thread_size(void);
int reel_thread_rank(void);
/* call this first thing in your thread_main */
void reel_thread_init(void* in);
void reel_thread_barrier(void);
void reel_thread_lock(void);
void reel_thread_unlock(void);

void reel_protect(void);

#endif
