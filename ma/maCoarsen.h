/****************************************************************************** 

  Copyright 2013 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#ifndef MA_COARSEN_H
#define MA_COARSEN_H

namespace ma {

class Adapt;

void coarsenOnce(Adapt* a);
bool coarsenMultiple(Adapt* a);
bool coarsen(Adapt* a);
bool coarsenLayer(Adapt* a);

void checkAllEdgeCollapses(Adapt* a, int modelDimension);
void findIndependentSet(Adapt* a);

int collapseAllEdges(Adapt* a, int modelDimension);

}

#endif
