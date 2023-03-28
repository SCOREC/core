/****************************************************************************** 

  Copyright 2011 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/
#ifndef PCU_PMPI_H
#define PCU_PMPI_H

#include "pcu_mpi.h"

#include <stdbool.h>

pcu_mpi_t* pcu_pmpi_init(MPI_Comm comm);
void pcu_pmpi_finalize(pcu_mpi_t** m);
int pcu_pmpi_size(pcu_mpi_t* self);
int pcu_pmpi_rank(pcu_mpi_t* self);
void pcu_pmpi_send(pcu_mpi_t*, pcu_message* m, MPI_Comm comm);
bool pcu_pmpi_receive(pcu_mpi_t*, pcu_message* m, MPI_Comm comm);
bool pcu_pmpi_done(pcu_mpi_t*, pcu_message* m);

#endif
