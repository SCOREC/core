/****************************************************************************** 

  Copyright 2011 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/
#include "pcu_pmpi.h"
#include "pcu_buffer.h"
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

void pcu_pmpi_send2(const pcu_mpi_t* self, pcu_message* m, int tag, MPI_Comm comm);
bool pcu_pmpi_receive2(const pcu_mpi_t*, pcu_message* m, int tag, MPI_Comm comm);

void pcu_pmpi_init(MPI_Comm comm, pcu_mpi_t* self)
{
  self->original_comm = comm;
  MPI_Comm_dup(comm,&(self->user_comm));
  MPI_Comm_dup(comm,&(self->coll_comm));
  MPI_Comm_size(comm,&(self->size));
  MPI_Comm_rank(comm,&(self->rank));
}

void pcu_pmpi_finalize(pcu_mpi_t* self)
{
  MPI_Comm_free(&(self->user_comm));
  MPI_Comm_free(&(self->coll_comm));
}

int pcu_pmpi_free(MPI_Comm* comm)
{
  return MPI_Comm_free(comm);
}

int pcu_pmpi_split(MPI_Comm comm, int color, int key, MPI_Comm* newcomm)
{
  return MPI_Comm_split(comm,color,key,newcomm);
}

int pcu_pmpi_size(const pcu_mpi_t* self)
{
  return self->size;
}

int pcu_pmpi_rank(const pcu_mpi_t* self)
{
  return self->rank;
}

void pcu_pmpi_send(const pcu_mpi_t* self, pcu_message* m, MPI_Comm comm)
{
  pcu_pmpi_send2(self, m,0,comm);
}

void pcu_pmpi_send2(const pcu_mpi_t* self, pcu_message* m, int tag, MPI_Comm comm)
{
  // silence warning
  (void)self;
  if( m->buffer.size > (size_t)INT_MAX ) {
    fprintf(stderr, "ERROR PCU message size exceeds INT_MAX... exiting\n");
    abort();
  }
  MPI_Issend(
      m->buffer.start,
      (int)(m->buffer.size),
      MPI_BYTE,
      m->peer,
      tag,
      comm,
      &(m->request));
}

bool pcu_pmpi_done(const pcu_mpi_t* self, pcu_message* m)
{
  // silence warning
  (void)self;
  int flag;
  MPI_Test(&(m->request),&flag,MPI_STATUS_IGNORE);
  return flag;
}

bool pcu_pmpi_receive(const pcu_mpi_t* self, pcu_message* m, MPI_Comm comm)
{
  // silence warning
  (void)self;
  return pcu_pmpi_receive2(self, m,0,comm);
}

bool pcu_pmpi_receive2(const pcu_mpi_t* self, pcu_message* m, int tag, MPI_Comm comm)
{
  // silence warning
  (void)self;
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

int  pcu_pmpi_allreduce(const void* sendbuf, void* recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
{
   return MPI_Allreduce(sendbuf,recvbuf,count, datatype, op, comm);
}

int  pcu_pmpi_allgather(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm)
{
  return MPI_Allgather(sendbuf,sendcount,sendtype,recvbuf,recvcount,recvtype,comm);
}
