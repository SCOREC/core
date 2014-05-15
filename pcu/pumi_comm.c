/****************************************************************************** 

  (c) 2011-2014 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/
#include <pumi_comm.h>
#include "PCU.h"

#define CALL(f) return (f)?PUMI_COMM_FAIL:0

int PUMI_Thrd_Run(int nthreads, PUMI_Thrd_Func function, void** in_out)
{
  CALL(PCU_Thrd_Run(nthreads,function,in_out));
}

int PUMI_Thrd_Rank(int* rank)
{
  *rank = PCU_Thrd_Self();
  return 0;
}

int PUMI_Thrd_Size(int* size)
{
  *size = PCU_Thrd_Peers();
  return 0;
}

int PUMI_Comm_Rank(int* rank)
{
  CALL(PCU_Comm_Rank(rank));
}

int PUMI_Comm_Size(int* size)
{
  CALL(PCU_Comm_Size(size));
}

int PUMI_Comm_Init(void)
{
  CALL(PCU_Comm_Init());
}

int PUMI_Comm_Start(int method)
{
  CALL(PCU_Comm_Start((PCU_Method)method));
}

int PUMI_Comm_Pack(int to_rank, const void* data, size_t size)
{
  CALL(PCU_Comm_Pack(to_rank,data,size));
}

int PUMI_Comm_Packed(int to_rank, size_t* size)
{
  CALL(PCU_Comm_Packed(to_rank,size));
}

int PUMI_Comm_Write(int to_rank, const void* data, size_t size)
{
  CALL(PCU_Comm_Write(to_rank,data,size));
}

int PUMI_Comm_Send(void)
{
  CALL(PCU_Comm_Send());
}

int PUMI_Comm_Receive(bool* done)
{
  CALL(PCU_Comm_Receive(done));
}

bool PUMI_Comm_Read(int* from_rank, void** data, size_t* size)
{
  return PCU_Comm_Read(from_rank,data,size);
}

int PUMI_Comm_From(int* from_rank)
{
  CALL(PCU_Comm_From(from_rank));
}

int PUMI_Comm_Received(size_t* size)
{
  CALL(PCU_Comm_Received(size));
}

int PUMI_Comm_Unpack(void* data, size_t size)
{
  CALL(PCU_Comm_Unpack(data,size));
}

bool PUMI_Comm_Unpacked(void)
{
  return PCU_Comm_Unpacked();
}

int PUMI_Comm_Free(void)
{
  CALL(PCU_Comm_Free());
}

