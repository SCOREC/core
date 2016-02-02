/****************************************************************************** 

  Copyright 2015 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/
#include "reel.h"
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <stddef.h>

void reel_fail(const char* format, ...)
{
  fprintf(stderr, "PUMI error: ");
  va_list ap;
  va_start(ap, format);
  vfprintf(stderr, format, ap);
  va_end(ap);
  fprintf(stderr, "\n");
  abort();
}

#if defined(__linux__) || defined(__APPLE__)
#include <execinfo.h> /* backtrace for pcu_trace */
#include <signal.h> /* signal for pcu_protect */

void reel_trace(void)
{
  static void* buf[64];
  int n;
  n = backtrace(buf, 64);
  backtrace_symbols_fd(buf, n, 2);
  fflush(stderr);
}

static void catch(int s)
{
  static volatile sig_atomic_t already_dying = 0;
  if (already_dying)
    return;
  already_dying = 1;
  fprintf(stderr, "signal %d caught by pcu\n", s);
  reel_trace();
  signal(s, SIG_DFL);
  raise(s);
}

void reel_protect(void)
{
  signal(SIGABRT, catch);
  signal(SIGSEGV, catch);
  signal(SIGINT, catch);
  signal(SIGFPE, catch);
}
#else
void reel_protect(void)
{
  reel_fail(stderr, "reel_protect only supported on Linux and OS X");
}
#endif
