#ifndef APF_AGGREGATE_NUMBERING_H_
#define APF_AGGREGATE_NUMBERING_H_
#include "apfNew.h"
#include "apfMesh.h"
#include "apfNumbering.h"
#include <PCU.h>
namespace apf
{
  // class declaration
  template <class T>
  class AggregateNumberingOf;
  typedef AggregateNumberingOf<int> AggNumbering;
  typedef AggregateNumberingOf<long> GlobalAggNumbering;
  // numbering creation functions
  AggNumbering * createAggNumbering(Field * f,
                                    int blocks,
                                    int dofs_per_block,
                                    MPI_Comm cm,
                                    Sharing * share = NULL);
  AggNumbering * createAggNumbering(Mesh * m,
                                    const char * name,
                                    FieldShape * shape,
                                    int blocks,
                                    int dofs_per_block,
                                    MPI_Comm cm,
                                    Sharing * share = NULL);
  /*
   * @brief Since the class implementation is not public
   *        we need an API to up-cast AggNumbering to
   *        Numbering.
   */
  Numbering * getNumbering(AggNumbering * n);
  int countBlocks(AggNumbering * n);
  int countDOFsPerBlock(AggNumbering * n);
}
#endif
/******************************************************************************
  Copyright 2018 Scientific Computation Research Center,
      Rensselaer Polytechnic Institute. All rights reserved.
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
*******************************************************************************/
