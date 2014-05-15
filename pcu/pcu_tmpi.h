/****************************************************************************** 

  (c) 2011-2014 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/
#ifndef PCU_TMPI_H
#define PCU_TMPI_H

#include "pcu_mpi.h"

int pcu_tmpi_size(void);
int pcu_tmpi_rank(void);
void pcu_tmpi_send(pcu_message* m, MPI_Comm comm);
bool pcu_tmpi_done(pcu_message* m);
bool pcu_tmpi_receive(pcu_message* m, MPI_Comm comm);

void pcu_tmpi_check_support(void);

extern pcu_mpi pcu_tmpi;

#endif
