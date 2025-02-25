/****************************************************************************** 

  Copyright 2011 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/
#ifndef PCU_MPI_H
#define PCU_MPI_H

#include "pcu_defines.h"
#include "pcu_buffer.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
  pcu_buffer buffer;
  PCU_Request request;
  int peer;
} pcu_message;

void pcu_make_message(pcu_message* m);
void pcu_free_message(pcu_message* m);

struct pcu_mpi_struct
{
  PCU_Comm user_comm;
  PCU_Comm coll_comm;
  int rank;
  int size;
};
typedef struct pcu_mpi_struct pcu_mpi_t;

int pcu_mpi_size(const pcu_mpi_t*);
int pcu_mpi_rank(const pcu_mpi_t*);
void pcu_mpi_send(const pcu_mpi_t*, pcu_message* m, PCU_Comm comm);
bool pcu_mpi_done(const pcu_mpi_t*, pcu_message* m);
bool pcu_mpi_receive(const pcu_mpi_t*, pcu_message* m, PCU_Comm comm);
void pcu_mpi_init(PCU_Comm comm, pcu_mpi_t* mpi);
void pcu_mpi_finalize(pcu_mpi_t* mpi);
int  pcu_mpi_split(const pcu_mpi_t* mpi, int color, int key, PCU_Comm* newcomm);
int  pcu_mpi_dup(const pcu_mpi_t* mpi, PCU_Comm* newcomm);

#ifdef __cplusplus
}
#endif
#endif
