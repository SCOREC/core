/****************************************************************************** 

  Copyright 2013 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#ifndef MA_H
#define MA_H

/** \page ma MeshAdapt
   This is the main API for SCOREC's MeshAdapt software, version 2

   These functions run parallel mesh adaptation based on various
   types of size fields, with optional user-defined solution transfer.

   If no solution transfer is specified, MeshAdapt will automatically
   transfer all nodal fields and some integration point fields.

   Features such a snapping and matched entity support are activated
   automatically if the input structures support them.

   An interface is also provided where users can pick and choose
   what MeshAdapt does and adjust some parameters.

 - The main API is in ma.h
 - Most users may need to also tweak settings using maInput.h
 - For defining size functions, see maSize.h
 - In you need special solution transfer, see maSolutionTransfer.h
*/

/** \file ma.h
  \brief The MeshAdapt interface */

#include "maInput.h"

/** \namespace ma
    \brief All MeshAdapt symbols */
namespace ma {

/** \brief adapt based on an isotropic function
   \details see maSize.h for how to define a function */
void adapt(Mesh* m, IsotropicFunction* f, SolutionTransfer* s=0);
/** \brief adapt based on an anisotropic function */
void adapt(Mesh* m, AnisotropicFunction* f, SolutionTransfer* s=0);
/** \brief adapt with custom configuration
  \details see maInput.h for details.
  note that this function will delete the Input object */
void adapt(Input* in);
/** \brief run uniform refinement, plus snapping and shape correction */
void runUniformRefinement(Mesh* m, int n=1, SolutionTransfer* s=0);
/** \brief run uniform refinement with matched entity support
  \details currently this supports snapping but not shape correction */
void adaptMatching(Mesh* m, int n=1, SolutionTransfer* s=0);

}

#endif
