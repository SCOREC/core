/****************************************************************************** 

  Copyright 2011 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/
#include "pcu_mpi.h"
#include <assert.h>

static pcu_mpi* global_mpi;

void pcu_make_message(pcu_message* m)
{
  pcu_make_buffer(&(m->buffer));
}

void pcu_free_message(pcu_message* m)
{
  pcu_free_buffer(&(m->buffer));
}

void pcu_set_mpi(pcu_mpi* m)
{
  global_mpi = m;
}

pcu_mpi* pcu_get_mpi(void)
{
  return global_mpi;
}

int pcu_mpi_size(void)
{
  return global_mpi->size();
}

int pcu_mpi_rank(void)
{
  return global_mpi->rank();
}

static void check_rank(int rank)
{
  assert(0 <= rank);
  assert(rank < pcu_mpi_size());
}

void pcu_mpi_send(pcu_message* m, MPI_Comm comm)
{
  check_rank(m->peer);
  global_mpi->send(m,comm);
}

bool pcu_mpi_done(pcu_message* m)
{
  return global_mpi->done(m);
}

bool pcu_mpi_receive(pcu_message* m, MPI_Comm comm)
{
  if (m->peer != MPI_ANY_SOURCE)
    check_rank(m->peer);
  return global_mpi->receive(m,comm);
}
