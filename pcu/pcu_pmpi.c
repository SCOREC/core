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
  self->owned = 0;
}

void pcu_pmpi_finalize(pcu_mpi_t* self)
{
  MPI_Comm_free(&(self->user_comm));
  MPI_Comm_free(&(self->coll_comm));
  // Prevent accidental freeing of MPI_COMM_WORLD.
  int result;
  MPI_Comm_compare(self->original_comm, MPI_COMM_WORLD, &result);
  if (self->owned && result != MPI_IDENT)
    MPI_Comm_free(&(self->original_comm));
}

int pcu_pmpi_split(const pcu_mpi_t *mpi, int color, int key, MPI_Comm* newcomm)
{
  return MPI_Comm_split(mpi->original_comm,color,key,newcomm);
}

int pcu_pmpi_dup(const pcu_mpi_t *mpi, PCU_Comm* newcomm)
{
  return MPI_Comm_dup(mpi->user_comm, newcomm);
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

