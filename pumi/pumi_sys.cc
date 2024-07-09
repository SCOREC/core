/****************************************************************************** 

  (c) 2004-2017 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
 
*******************************************************************************/
#include "pumi.h"
#include <lionPrint.h>
#include <mpi.h>

//************************************
//************************************
//      0- SYSTEM-LEVEL FUNCTIONS
//************************************
//************************************
void pumi_load_pcu(pcu::PCU *PCUObj){
  pumi::instance()->initializePCU(PCUObj);
}

int pumi_size()
{
  return pumi::instance()->getPCU()->Peers();
}

int pumi_rank()
{
  return pumi::instance()->getPCU()->Self();
}

void pumi_sync()
{
  MPI_Barrier(pumi::instance()->getPCU()->GetMPIComm());
}

#include <sys/utsname.h>
#include <sys/resource.h>
void pumi_printSys()
{
  if (pumi::instance()->getPCU()->Self()) return;
  struct utsname u;
  if (uname(&u) == 0)
    lion_oprint(1,"[%s] %s %s %s %s %s\n\n",
           __func__, u.sysname, u.nodename, u.release, u.version, u.machine);
  fflush(stdout);
}

double pumi_getTime()
{
  struct rusage ruse_now;
  getrusage(RUSAGE_SELF, &ruse_now);
  return double(ruse_now.ru_utime.tv_sec) + double(ruse_now.ru_utime.tv_usec)/1000000.0;
}

double pumi_getMem()
{
  return pcu::GetMem();
}

void pumi_printTimeMem(const char* msg, double time, double memory)
{
  if (!pumi::instance()->getPCU()->Self())
  {
    lion_oprint(1,"%-20s %6.3f sec %7.3f MB \n", msg, time, memory);
    fflush(stdout);
  }
}

