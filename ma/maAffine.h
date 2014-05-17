/****************************************************************************** 

  Copyright (c) 2004-2014 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/

#ifndef MA_AFFINE_H
#define MA_AFFINE_H

#include "apfMatrix.h"

namespace ma {

class Affine
{
  public:
    apf::Vector3 operator*(apf::Vector3 const& x) const
    {
      return A * x + b;
    }
    apf::Matrix3x3 A; 
    apf::Vector3 b;
};

inline Affine invert(Affine const& a)
{
  Affine ainv;
  ainv.A = apf::invert(a.A);
  ainv.b = (ainv.A * a.b) * -1;
  return ainv;
}

}

#endif
