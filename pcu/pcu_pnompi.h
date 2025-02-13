/****************************************************************************** 

  Copyright 2011 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/
#ifndef PCU_NO_PMPI_H
#define PCU_NO_PMPI_H

#include "pcu_mpi.h"

#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

void pcu_pmpi_init(MPI_Comm comm, pcu_mpi_t *mpi);
void pcu_pmpi_finalize(pcu_mpi_t *m);
int pcu_pmpi_size(const pcu_mpi_t *self);
int pcu_pmpi_rank(const pcu_mpi_t *self);
void pcu_pmpi_send(const pcu_mpi_t*, pcu_message *m, MPI_Comm comm);
bool pcu_pmpi_receive(const pcu_mpi_t*, pcu_message *m, MPI_Comm comm);
bool pcu_pmpi_done(const pcu_mpi_t*, pcu_message *m);

int  pcu_pmpi_split(MPI_Comm c, int color, int key, MPI_Comm* nc);
int  pcu_pmpi_free(MPI_Comm* comm);
int  pcu_pmpi_barrier(MPI_Comm comm);

double MPI_Wtime(void);

#ifdef __cplusplus
}
#endif

#endif
