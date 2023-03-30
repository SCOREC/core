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

#define PCU_DO_PRAGMA_(x) _Pragma (#x)
#define PCU_DO_PRAGMA(x) PCU_DO_PRAGMA_(x)

#define PCU_IGNORE_DIAGNOSTIC_START(x) \
  PCU_DO_PRAGMA(GCC diagnostic push) \
  PCU_DO_PRAGMA(GCC diagnostic ignored #x) \
  PCU_DO_PRAGMA(clang diagnostic push) \
  PCU_DO_PRAGMA(clang diagnostic ignored #x)

#define PCU_IGNORE_DIAGNOSTIC_END \
  PCU_DO_PRAGMA(GCC diagnostic pop) \
  PCU_DO_PRAGMA(clang diagnostic pop)

#define PCU_ALWAYS_ASSERT(cond) \
  PCU_IGNORE_DIAGNOSTIC_START(-Wdeprecated-declarations)      \
  do {                                                      \
    if (! (cond)) {                                         \
      char omsg[2048];                                      \
      sprintf(omsg, "%s failed at %s + %d \n",              \
        #cond, __FILE__, __LINE__);                         \
      PCU_Assert_Fail(omsg);                                \
    }                                                       \
  } while (0)                                               \
  PCU_IGNORE_DIAGNOSTIC_END
#define PCU_ALWAYS_ASSERT_VERBOSE(cond, msg)      \
  PCU_IGNORE_DIAGNOSTIC_START(-Wdeprecated-declarations)      \
  do {                                            \
    if (! (cond)) {                               \
      char omsg[2048];                            \
      sprintf(omsg, "%s failed at %s + %d \n %s", \
        #cond, __FILE__, __LINE__, msg);          \
      PCU_Assert_Fail(omsg);                      \
    }                                             \
  } while(0)                                      \
  PCU_IGNORE_DIAGNOSTIC_END

#ifdef NDEBUG
#define PCU_DEBUG_ASSERT(cond)
#define PCU_DEBUG_ASSERT_VERBOSE(cond, msg)
#else
#define PCU_DEBUG_ASSERT(cond) \
  PCU_ALWAYS_ASSERT(cond)
#define PCU_DEBUG_ASSERT_VERBOSE(cond, msg) \
  PCU_ALWAYS_ASSERT_VERBOSE(cond, msg)
#endif

#endif
