/****************************************************************************** 

  Copyright 2013 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#ifndef MA_H
#define MA_H

#include "maInput.h"

namespace ma {

/* This is the main API for SCOREC's MeshAdapt software, version 2

   These functions run parallel mesh adaptation based on various
   types of size fields, with optional user-specific solution transfer.

   If no solution transfer is specified, MeshAdapt will automatically
   transfer all linear nodal fields.

   Features such a snapping and matched entity support are activated
   automatically if the input structures support them.

   An interface is also provided where users can pick and choose
   what MeshAdapt does and adjust some parameters.
*/

/* adapt based on an isotropic function
   see maSize.h for how to define a function */
void adapt(Mesh* m, IsotropicFunction* f, SolutionTransfer* s=0);
/* adapt based on an anisotropic function */
void adapt(Mesh* m, AnisotropicFunction* f, SolutionTransfer* s=0);
/* adapt with custom configuration
   see maInput.h for details
   note that this function will delete the Input object */
void adapt(Input* in);
/* run uniform refinement, plus snapping and shape correction if applicable */
void runUniformRefinement(Mesh* m, int n=1, SolutionTransfer* s=0);
/* run uniform refinement with matched entity support
   currently this supports snapping but not shape correction */
void adaptMatching(Mesh* m, int n=1, SolutionTransfer* s=0);

}

#endif
