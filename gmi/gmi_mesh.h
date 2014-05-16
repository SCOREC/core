/****************************************************************************** 

  (c) 2011-2014 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/
#ifndef GMI_MESH_H
#define GMI_MESH_H

#include "gmi.h"

void gmi_register_mesh(void);

struct gmi_mesh {
  struct gmi_model model;
  int* tags[4];
};

struct gmi_iter* gmi_mesh_begin(struct gmi_model* m, int dim);

#endif

