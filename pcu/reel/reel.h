/****************************************************************************** 

  Copyright 2015 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/
#ifndef REEL_H
#define REEL_H

#ifdef __cplusplus
extern "C" {
#endif

void reel_fail(const char* format, ...)
  __attribute__((noreturn,format(printf,1,2)));

void reel_protect(void);

#ifdef __cplusplus
}
#endif

#endif
