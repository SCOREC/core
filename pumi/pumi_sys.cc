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

void pumi_finalize(bool do_mpi_finalize)
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
  MPI_Barrier(MPI_COMM_WORLD);
}

#include <sys/utsname.h>
#include <sys/resource.h>
void pumi_info()
{
  if (PCU_Comm_Self()) return;
  struct utsname u;
  if (uname(&u) == 0)
    printf("[%s] %s %s %s %s %s\n\n",
           __func__, u.sysname, u.nodename, u.release, u.version, u.machine);
  fflush(stdout);
}
