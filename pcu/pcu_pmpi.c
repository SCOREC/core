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

void pcu_pmpi_send2(pcu_mpi_t* self, pcu_message* m, int tag, MPI_Comm comm);
bool pcu_pmpi_receive2(pcu_mpi_t*, pcu_message* m, int tag, MPI_Comm comm);
void pcu_pmpi_destruct(pcu_mpi_t* self);
void pcu_pmpi_construct(pcu_mpi_t* self, MPI_Comm comm);


pcu_mpi_t* pcu_pmpi_init(MPI_Comm comm)
{
  pcu_mpi_t* m = (pcu_mpi_t*)malloc(sizeof(pcu_mpi_t));
  pcu_pmpi_construct(m,comm);
  return m;
}

void pcu_pmpi_finalize(pcu_mpi_t** m)
{
  pcu_pmpi_destruct(*m);
  free(*m);
  *m = NULL;
}
void pcu_pmpi_destruct(pcu_mpi_t* self) {
  MPI_Comm_free(&(self->user_comm));
  MPI_Comm_free(&(self->coll_comm));
}
void pcu_pmpi_construct(pcu_mpi_t* self, MPI_Comm comm) {
  self->original_comm = comm;
  MPI_Comm_dup(comm,&(self->user_comm));
  MPI_Comm_dup(comm,&(self->coll_comm));
  MPI_Comm_size(comm,&(self->size));
  MPI_Comm_rank(comm,&(self->rank));
}

int pcu_pmpi_size(pcu_mpi_t* self)
{
  return self->size;
}

int pcu_pmpi_rank(pcu_mpi_t* self)
{
  return self->rank;
}

void pcu_pmpi_send(pcu_mpi_t* self, pcu_message* m, MPI_Comm comm)
{
  pcu_pmpi_send2(self, m,0,comm);
}

void pcu_pmpi_send2(pcu_mpi_t* self, pcu_message* m, int tag, MPI_Comm comm)
{
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

bool pcu_pmpi_done(pcu_mpi_t* self, pcu_message* m)
{
  int flag;
  MPI_Test(&(m->request),&flag,MPI_STATUS_IGNORE);
  return flag;
}

bool pcu_pmpi_receive(pcu_mpi_t* self, pcu_message* m, MPI_Comm comm)
{
  return pcu_pmpi_receive2(self, m,0,comm);
}

bool pcu_pmpi_receive2(pcu_mpi_t* self, pcu_message* m, int tag, MPI_Comm comm)
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
