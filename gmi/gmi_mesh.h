/****************************************************************************** 

  Copyright 2014 Scientific Computation Research Center, 
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

void gmi_write_dmg(struct gmi_model* m, const char* filename);

struct gmi_mesh {
  struct gmi_model model;
  int* tags[4];
  struct gmi_set** up[4];
  struct gmi_set** down[4];
};

struct gmi_iter* gmi_mesh_begin(struct gmi_model* m, int dim);
struct gmi_ent* gmi_mesh_next(struct gmi_model* m, struct gmi_iter* it);
void gmi_mesh_end(struct gmi_model* m, struct gmi_iter* i);
int gmi_mesh_dim(struct gmi_model* m, struct gmi_ent* e);
int gmi_mesh_tag(struct gmi_model* m, struct gmi_ent* e);
void gmi_mesh_destroy(struct gmi_model* m);
int gmi_index(struct gmi_ent* e);
struct gmi_ent* gmi_identify(int dim, int idx);

#ifdef __cplusplus
}
#endif

#endif

