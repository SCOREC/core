/****************************************************************************** 

  Copyright 2014 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/
#ifndef GMI_BASE_H
#define GMI_BASE_H

#include "gmi.h"

#ifdef __cplusplus
extern "C" {
#endif

/* base struct for all the internal gmi structures:
   mesh, null, analytic, etc. */

struct gmi_base {
  struct gmi_model model;
  int* tags[4];
};

struct gmi_iter* gmi_base_begin(struct gmi_model* m, int dim);
struct gmi_ent* gmi_base_next(struct gmi_model* m, struct gmi_iter* it);
void gmi_base_end(struct gmi_model* m, struct gmi_iter* i);
int gmi_base_dim(struct gmi_model* m, struct gmi_ent* e);
int gmi_base_tag(struct gmi_model* m, struct gmi_ent* e);
struct gmi_ent* gmi_base_find(struct gmi_model* m, int dim, int tag);
void gmi_base_destroy(struct gmi_model* m);
int gmi_base_index(struct gmi_ent* e);
struct gmi_ent* gmi_base_identify(int dim, int idx);
void gmi_base_read_dmg(struct gmi_base* m, FILE* f);

#ifdef __cplusplus
}
#endif

#endif

