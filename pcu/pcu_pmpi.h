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

void pcu_pmpi_init(MPI_Comm comm);
void pcu_pmpi_finalize(void);
int  pcu_pmpi_size(void);
int  pcu_pmpi_rank(void);
void pcu_pmpi_send(pcu_message* m, MPI_Comm comm);
bool pcu_pmpi_receive(pcu_message* m, MPI_Comm comm);
void pcu_pmpi_send2(pcu_message* m, int tag, MPI_Comm comm);
bool pcu_pmpi_receive2(pcu_message* m, int tag, MPI_Comm comm);
bool pcu_pmpi_done(pcu_message* m);

void pcu_pmpi_switch(MPI_Comm new_comm);
MPI_Comm pcu_pmpi_comm(void);
int  pcu_pmpi_lcl_to_frn(int lcl, MPI_Comm frn_cm);
int  pcu_pmpi_frn_to_lcl(int frn, MPI_Comm frn_cm);

extern pcu_mpi pcu_pmpi;

extern MPI_Comm pcu_user_comm;
extern MPI_Comm pcu_coll_comm;

#endif
