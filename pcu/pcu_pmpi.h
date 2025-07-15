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
#ifdef __cplusplus
extern "C" {
#endif
void pcu_pmpi_init(PCU_Comm comm, pcu_mpi_t *mpi);
void pcu_pmpi_finalize(pcu_mpi_t *m);
int pcu_pmpi_size(const pcu_mpi_t *self);
int pcu_pmpi_rank(const pcu_mpi_t *self);
void pcu_pmpi_send(const pcu_mpi_t *, pcu_message *m, PCU_Comm comm);
bool pcu_pmpi_receive(const pcu_mpi_t *, pcu_message *m, PCU_Comm comm);
bool pcu_pmpi_done(const pcu_mpi_t *, pcu_message *m);

int pcu_pmpi_split(const pcu_mpi_t *, int color, int key, PCU_Comm* newcomm);
int pcu_pmpi_dup(const pcu_mpi_t *, PCU_Comm* newcomm);

#ifdef __cplusplus
}
#endif
#endif
