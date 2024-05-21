/****************************************************************************** 

  (c) 2004-2017 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
 
*******************************************************************************/
#include "pumi.h"
#include <lionPrint.h>

//************************************
//************************************
//      0- SYSTEM-LEVEL FUNCTIONS
//************************************
//************************************

void pumi_finalize(pcu::PCU *PCUObj, bool)
{
  delete PCUObj;
}

int pumi_size(pcu::PCU *PCUObj)
{
  return PCUObj->Peers();
}

int pumi_rank(pcu::PCU *PCUObj)
{
  return PCUObj->Self();
}

void pumi_sync(pcu::PCU *PCUObj)
{
  MPI_Barrier(PCUObj->GetMPIComm());
}

#include <sys/utsname.h>
#include <sys/resource.h>
void pumi_printSys(pcu::PCU *PCUObj)
{
  if (PCUObj->Self()) return;
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

void pumi_printTimeMem(const char* msg, double time, double memory, pcu::PCU *PCUObj)
{
  if (!PCUObj->Self())
  {
    lion_oprint(1,"%-20s %6.3f sec %7.3f MB \n", msg, time, memory);
    fflush(stdout);
  }
}

