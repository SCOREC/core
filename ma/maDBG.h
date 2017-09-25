/****************************************************************************** 

  Copyright 2013 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#ifndef MA_DBG_H
#define MA_DBG_H

#include <apfMesh2.h>
#include <apfShape.h>
#include <apfNumbering.h>
#include "maAdapt.h"

#include <vector>
#include <assert.h>

namespace ma_dbg {

void colorEntitiesOfDimWithQual(ma::Adapt* a, int dim, const std::vector<double> & quals, const char* fieldName);
void writeMesh(ma::Mesh* m, const char* prefix, const char* suffix);
void measureLinearQualties(ma::Adapt* a,
     std::vector<double> &lq, bool inMetric = true);
void dumpMeshWithQualities(ma::Adapt* a, int iter, const char* prefix);

}
#endif
