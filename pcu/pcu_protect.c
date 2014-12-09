#include <stdio.h>
#include "pcu_io.h" 

#if defined(__linux__) || defined(__APPLE__)
#include <execinfo.h> /* backtrace for pcu_trace */
#include <signal.h> /* signal for pcu_protect */

void PCU_Trace(void)
{
  FILE* fp = pcu_open_parallel("trace","txt");
  int fd = fileno(fp);
  static void* buf[64];
  int n;
  n = backtrace(buf, 64);
  backtrace_symbols_fd(buf, n, fd);
  fclose(fp);
}

static void catch(int s)
{
  static volatile sig_atomic_t already_dying = 0;
  if (already_dying)
    return;
  already_dying = 1;
  fprintf(stderr, "signal %d caught by pcu\n",s);
  PCU_Trace();
  signal(s, SIG_DFL);
  raise(s);
}

void PCU_Protect(void)
{
  signal(SIGABRT, catch);
  signal(SIGSEGV, catch);
  signal(SIGINT, catch);
}
#else
void PCU_Protect(void)
{
  fprintf(stderr,"PCU_Protect only supported on Linux and OS X\n");
}
#endif

