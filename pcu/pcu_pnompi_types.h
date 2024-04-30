/******************************************************************************

  Copyright 2011 Scientific Computation Research Center,
      Rensselaer Polytechnic Institute. All rights reserved.

  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/
#ifndef PCU_PNOMPI_TYPES_H
#define PCU_PNOMPI_TYPES_H

typedef int MPI_Comm;
typedef int MPI_Request;
typedef int MPI_Datatype;
typedef int MPI_Op;
#define MPI_COMM_WORLD 0
#define MPI_ANY_SOURCE -1
#define MPI_SUM 1
#define MPI_INT 2

#endif // PCU_PNOMPI_TYPES_H
