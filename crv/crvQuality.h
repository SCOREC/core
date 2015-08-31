/******************************************************************************

  Copyright 2015 Scientific Computation Research Center,
      Rensselaer Polytechnic Institute. All rights reserved.

  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.

*******************************************************************************/

#ifndef CRVQUALITY_H
#define CRVQUALITY_H

#include "crv.h"

namespace crv {

void subdivideBezierEdgeJacobianDet(int P,
    apf::NewArray<double>& nodes,
    apf::NewArray<double> (&subNodes)[2]);

void subdivideBezierTriangleJacobianDet(int P,
    apf::NewArray<double>& nodes,
    apf::NewArray<double> (&subNodes)[4]);

}

#endif
