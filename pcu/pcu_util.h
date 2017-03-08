/******************************************************************************

  Copyright 2011 Scientific Computation Research Center,
      Rensselaer Polytechnic Institute. All rights reserved.

  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/
#ifndef PCU_UTIL_H
#define PCU_UTIL_H

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

void PCU_Assert_Fail(const char* msg) __attribute__ ((noreturn));

#ifdef __cplusplus
} /* extern "C" */
#endif

#define PCU_ALWAYS_ASSERT(...)                    \
  do {                                            \
    if (! (__VA_ARGS__)) {                        \
      char omsg[512];                             \
      sprintf(omsg, "%s failed at %s + %d \n",    \
        #__VA_ARGS__, __FILE__, __LINE__);        \
      PCU_Assert_Fail(omsg);                      \
    }                                             \
  } while (0)

#define PCU_ALWAYS_ASSERT_VERBOSE(msg, ...)       \
  do {                                            \
    if (! (__VA_ARGS__)) {                        \
      char omsg[512];                             \
      sprintf(omsg, "%s failed at %s + %d \n %s", \
        #__VA_ARGS__, __FILE__, __LINE__, msg);   \
      PCU_Assert_Fail(omsg);                      \
    }                                             \
  } while(0)

#ifdef NDEBUG
#define PCU_DEBUG_ASSERT(...)
#define PCU_DEBUG_ASSERT_VERBOSE(msg, ...)
#else
#define PCU_DEBUG_ASSERT(...) \
  PCU_ALWAYS_ASSERT(__VA_ARGS__)
#define PCU_DEBUG_ASSERT_VERBOSE(msg, ...) \
  PCU_ALWAYS_ASSERT_VERBOSE(msg, __VA_ARGS__)
#endif

#endif
