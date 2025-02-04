/****************************************************************************** 

  Copyright 2011 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/
#include "pcu_mpi.h"
#include "pcu_util.h"
#include "pcu_pmpi.h"

void pcu_make_message(pcu_message* m)
{
  pcu_make_buffer(&(m->buffer));
}

void pcu_free_message(pcu_message* m)
{
  pcu_free_buffer(&(m->buffer));
}

int pcu_mpi_size(const pcu_mpi_t* self)
{
  return pcu_pmpi_size(self);
}

int pcu_mpi_rank(const pcu_mpi_t* self)
{
  return pcu_pmpi_rank(self);
}

static void check_rank(const pcu_mpi_t* self, int rank)
{
  (void)rank;
  PCU_ALWAYS_ASSERT(0 <= rank);
  PCU_ALWAYS_ASSERT(rank < pcu_mpi_size(self));
}

void pcu_mpi_send(const pcu_mpi_t* self, pcu_message* m, MPI_Comm comm)
{
  check_rank(self, m->peer);
  PCU_ALWAYS_ASSERT(comm == self->user_comm || comm == self->coll_comm);
  pcu_pmpi_send(self, m, comm);
}

bool pcu_mpi_done(const pcu_mpi_t* self, pcu_message* m)
{
  return pcu_pmpi_done(self, m);
}

bool pcu_mpi_receive(const pcu_mpi_t* self, pcu_message* m, MPI_Comm comm)
{
  if (m->peer != MPI_ANY_SOURCE)
    check_rank(self, m->peer);
  return pcu_pmpi_receive(self, m, comm);
}
void pcu_mpi_init(MPI_Comm comm, pcu_mpi_t* mpi) {
  pcu_pmpi_init(comm, mpi);
}
void pcu_mpi_finalize(pcu_mpi_t* mpi) {
  pcu_pmpi_finalize(mpi);
}
