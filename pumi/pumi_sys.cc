/****************************************************************************** 

  (c) 2004-2016 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
 
*******************************************************************************/
#include <PCU.h>
#include "pumi.h"
#include <mpi.h>

//************************************
//************************************
//      0- SYSTEM-LEVEL FUNCTIONS
//************************************
//************************************

void pumi_start()
{
  PCU_Comm_Init();
}

void pumi_finalize(bool)
{
  PCU_Comm_Free();
}

int pumi_size()
{
  return PCU_Comm_Peers();
}

int pumi_rank()
{
  return PCU_Comm_Self();
}

void pumi_sync(void)
{
  MPI_Barrier(PCU_Get_Comm());
}

#include <sys/utsname.h>
#include <sys/resource.h>
void pumi_printSys()
{
  if (PCU_Comm_Self()) return;
  struct utsname u;
  if (uname(&u) == 0)
    printf("[%s] %s %s %s %s %s\n\n",
           __func__, u.sysname, u.nodename, u.release, u.version, u.machine);
  fflush(stdout);
}

double pumi_getTime()
{
  struct rusage ruse_now;
  getrusage(RUSAGE_SELF, &ruse_now);
  return double(ruse_now.ru_utime.tv_sec) + double(ruse_now.ru_utime.tv_usec)/1000000.0;
}

#if defined(__APPLE__)
#include <mach/task.h>
#include <mach/mach_init.h>
#elif defined(__bgq__)
#include <spi/include/kernel/memory.h>
#else
#include <malloc.h> //warning - this is GNU-specific
#endif

double pumi_getMem()
{
  const double M = 1024*1024;
#if defined(__APPLE__)
  bool resident = true;
  struct task_basic_info t_info;
  mach_msg_type_number_t t_info_count = TASK_BASIC_INFO_COUNT;
  task_info(current_task(), TASK_BASIC_INFO, (task_info_t)&t_info, &t_info_count);
  size_t size = (resident ? t_info.resident_size : t_info.virtual_size);
  return (double)size/M;
#elif defined(__bgq__)
  size_t heap;
  Kernel_GetMemorySize(KERNEL_MEMSIZE_HEAP, &heap);
  return (double)heap/M;
#else
  struct mallinfo meminfo_now = mallinfo();
  return ((double)meminfo_now.arena)/M;
#endif
}

void pumi_printTimeMem(const char* msg, double time, double memory)
{
  if (!PCU_Comm_Self())
  {
    printf("%-20s %6.3f sec %7.3f MB \n", msg, time, memory);
    fflush(stdout);
  }
}

