/****************************************************************************** 

  Copyright 2014 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/
#ifndef GMI_BASE_H
#define GMI_BASE_H

#include "gmi.h"
#include "agm.h"

#ifdef __cplusplus
extern "C" {
#endif

struct gmi_lookup;

/* base struct for all the internal gmi structures:
   mesh, null, analytic, etc. */

struct gmi_base {
  struct gmi_model model;
  struct agm* topo;
  struct gmi_lookup* lookup;
};

struct gmi_ent* gmi_from_agm(struct agm_ent e);
struct agm_ent agm_from_gmi(struct gmi_ent* e);
struct agm* gmi_base_topo(struct gmi_model* m);

void gmi_base_init(struct gmi_base* m);
void gmi_base_reserve(struct gmi_base* m, int dim, int n);
void gmi_base_destroy(struct gmi_model* m);
struct gmi_iter* gmi_base_begin(struct gmi_model* m, int dim);
struct gmi_ent* gmi_base_next(struct gmi_model* m, struct gmi_iter* it);
void gmi_base_end(struct gmi_model* m, struct gmi_iter* i);
int gmi_base_dim(struct gmi_model* m, struct gmi_ent* e);
int gmi_base_tag(struct gmi_model* m, struct gmi_ent* e);
struct gmi_ent* gmi_base_find(struct gmi_model* m, int dim, int tag);
struct gmi_set* gmi_base_adjacent(struct gmi_model* m, struct gmi_ent* e,
    int dim);

void gmi_base_freeze(struct gmi_model* m);
void gmi_base_unfreeze(struct gmi_model* m);

int gmi_base_index(struct gmi_ent* e);
struct gmi_ent* gmi_base_identify(int dim, int idx);

void gmi_base_read_dmg(struct gmi_base* m, FILE* f);
void gmi_base_read_tess(struct gmi_base* m, FILE* f);

void gmi_base_set_tag(struct gmi_model* m, struct gmi_ent* e, int tag);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif

