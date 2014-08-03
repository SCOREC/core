/****************************************************************************** 

  Copyright 2011 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/
#include "pcu_pmpi.h"
#include "pcu_common.h"
#include "pcu_memory.h"

static int global_size;
static int global_rank;

MPI_Comm original_comm;
MPI_Comm pcu_user_comm;
MPI_Comm pcu_coll_comm;

pcu_mpi pcu_pmpi =
{ .size = pcu_pmpi_size,
  .rank = pcu_pmpi_rank,
  .send = pcu_pmpi_send,
  .done = pcu_pmpi_done,
  .receive = pcu_pmpi_receive };

void pcu_pmpi_init(MPI_Comm comm)
{
  original_comm = comm;
  MPI_Comm_dup(comm,&pcu_user_comm);
  MPI_Comm_dup(comm,&pcu_coll_comm);
  MPI_Comm_size(comm,&global_size);
  MPI_Comm_rank(comm,&global_rank);
}

void pcu_pmpi_finalize(void)
{
  MPI_Comm_free(&pcu_user_comm);
  MPI_Comm_free(&pcu_coll_comm);
}

int pcu_pmpi_size(void)
{
  return global_size;
}

int pcu_pmpi_rank(void)
{
  return global_rank;
}

void pcu_pmpi_send(pcu_message* m, MPI_Comm comm)
{
  pcu_pmpi_send2(m,0,comm);
}

void pcu_pmpi_send2(pcu_message* m, int tag, MPI_Comm comm)
{
  MPI_Issend(
      m->buffer.start,
      (int)(m->buffer.size),
      MPI_BYTE,
      m->peer,
      tag,
      comm,
      &(m->request));
}

bool pcu_pmpi_done(pcu_message* m)
{
  int flag;
  MPI_Test(&(m->request),&flag,MPI_STATUS_IGNORE);
  return flag;
}

bool pcu_pmpi_receive(pcu_message* m, MPI_Comm comm)
{
  return pcu_pmpi_receive2(m,0,comm);
}

bool pcu_pmpi_receive2(pcu_message* m, int tag, MPI_Comm comm)
{
  MPI_Status status;
  int flag;
  MPI_Iprobe(m->peer,tag,comm,&flag,&status);
  if (!flag)
    return false;
  m->peer = status.MPI_SOURCE;
  int count;
  MPI_Get_count(&status,MPI_BYTE,&count);
  pcu_resize_buffer(&(m->buffer),(size_t)count);
  MPI_Recv(
      m->buffer.start,
      (int)(m->buffer.size),
      MPI_BYTE,
      m->peer,
      tag,
      comm,
      MPI_STATUS_IGNORE);
  return true;
}

void pcu_pmpi_switch(MPI_Comm new_comm)
{
  pcu_pmpi_finalize();
  pcu_pmpi_init(new_comm);
}

MPI_Comm pcu_pmpi_comm(void)
{
  return original_comm;
}

