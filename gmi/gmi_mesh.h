/****************************************************************************** 

  Copyright 2014 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/
#ifndef GMI_MESH_H
#define GMI_MESH_H

/** \file gmi_mesh.h
  \brief GMI meshmodel interface */

#include "gmi_base.h"

#ifdef __cplusplus
extern "C" {
#endif

/** \brief register the meshmodel reader for .dmg files */
void gmi_register_mesh(void);

#ifdef __cplusplus
}
#endif

#endif

