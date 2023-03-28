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

int pcu_mpi_size(pcu_mpi_t* self)
{
  return pcu_pmpi_size(self);
}

int pcu_mpi_rank(pcu_mpi_t* self)
{
  return pcu_pmpi_rank(self);
}

static void check_rank(pcu_mpi_t* self, int rank)
{
  (void)rank;
  PCU_ALWAYS_ASSERT(0 <= rank);
  PCU_ALWAYS_ASSERT(rank < pcu_mpi_size(self));
}

void pcu_mpi_send(pcu_mpi_t* self, pcu_message* m, MPI_Comm comm)
{
  check_rank(self, m->peer);
  // verify that the communicator is one of the user or collective communicators
  //int user_result, coll_result;
  //MPI_Comm_compare(comm, self->user_comm, &user_result);
  //MPI_Comm_compare(comm, self->coll_comm, &coll_result);
  //PCU_ALWAYS_ASSERT(user_result == MPI_IDENT || coll_result == MPI_IDENT);
  PCU_ALWAYS_ASSERT(comm == self->user_comm || comm == self->coll_comm);
  pcu_pmpi_send(self, m, comm);
}

bool pcu_mpi_done(pcu_mpi_t* self, pcu_message* m)
{
  return pcu_pmpi_done(self, m);
}

bool pcu_mpi_receive(pcu_mpi_t* self, pcu_message* m, MPI_Comm comm)
{
  if (m->peer != MPI_ANY_SOURCE)
    check_rank(self, m->peer);
  return pcu_pmpi_receive(self, m, comm);
}
pcu_mpi_t* pcu_mpi_init(MPI_Comm comm) {
  return pcu_pmpi_init(comm);
}
void pcu_mpi_finalize(pcu_mpi_t** mpi) {
  pcu_pmpi_finalize(mpi);
}
