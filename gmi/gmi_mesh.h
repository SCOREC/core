/****************************************************************************** 

  (c) 2011-2014 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/
#ifndef GMI_MESH_H
#define GMI_MESH_H

#include "gmi.h"

#ifdef __cplusplus
extern "C" {
#endif

void gmi_register_mesh(void);

struct gmi_mesh {
  struct gmi_model model;
  int* tags[4];
};

#ifdef __cplusplus
}
#endif

#endif

