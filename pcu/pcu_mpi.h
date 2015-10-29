/****************************************************************************** 

  Copyright 2011 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/
#ifndef PCU_MPI_H
#define PCU_MPI_H

#include "pcu_buffer.h"
#include <mpi.h>

typedef struct
{
  pcu_buffer buffer;
  MPI_Request request;
  int peer;
} pcu_message;

void pcu_make_message(pcu_message* m);
void pcu_free_message(pcu_message* m);

typedef struct
{
  int (*size)(void);
  int (*rank)(void);
  void (*send)(pcu_message* m, MPI_Comm comm);
  bool (*done)(pcu_message* m);
  bool (*receive)(pcu_message* m, MPI_Comm comm);
} pcu_mpi;

void pcu_set_mpi(pcu_mpi* m);
pcu_mpi* pcu_get_mpi(void);
int pcu_mpi_size(void);
int pcu_mpi_rank(void);
void pcu_mpi_send(pcu_message* m, MPI_Comm comm);
bool pcu_mpi_done(pcu_message* m);
bool pcu_mpi_receive(pcu_message* m, MPI_Comm comm);

#endif
