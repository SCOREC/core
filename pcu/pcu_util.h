/******************************************************************************

  Copyright 2011 Scientific Computation Research Center,
      Rensselaer Polytechnic Institute. All rights reserved.

  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/
#ifndef PCU_UTIL_H
#define PCU_UTIL_H

#ifdef __cplusplus
extern "C" {
#endif

void PCU_Assert_Fail(const char* msg) __attribute__ ((noreturn));

#ifdef __cplusplus
} /* extern "C" */
#endif

#define PCU_ASSERT_IMPL(cond, msg, ...)           \
  do {                                            \
    if (!(cond)) {                                \
      char omsg[512];                             \
      sprintf(omsg, "%s failed at %s + %d \n %s", \
        #cond, __FILE__, __LINE__, msg);          \
      PCU_Assert_Fail(omsg);                      \
    }                                             \
  } while (0)

#define PCU_ASSERT(...) PCU_ASSERT_IMPL(__VA_ARGS__, "")

#ifdef NDEBUG
#define PCU_EXPECT(...)
#else
#define PCU_EXPECT(...) PCU_ASSERT(__VA_ARGS__)
#endif

#define PCU_ALWAYS_ASSERT(cond) PCU_ASSERT(cond)
#define PCU_ALWAYS_ASSERT_VERBOSE(cond, msg) PCU_ASSERT(cond, msg)
#define PCU_DEBUG_ASSERT(cond) PCU_EXPECT(cond)
#define PCU_DEBUG_ASSERT_VERBOSE(cond, msg) PCU_EXPECT(cond, msg)

#endif
