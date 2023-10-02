/****************************************************************************** 

  Copyright 2014 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/
#include "gmi_null.h"
#include "gmi_base.h"
#include "gmi_lookup.h"
#include <stdlib.h>

struct gmi_ent* gmi_null_find(struct gmi_model* m, int dim, int tag)
{
  struct gmi_base* b;
  struct agm_ent e;
  enum agm_ent_type t;
  b = (struct gmi_base*)m;
  t = agm_type_from_dim(dim);
  e = gmi_look_up(b->lookup, t, tag);
  if (agm_ent_null(e)) {
    e = agm_add_ent(b->topo, t);
    gmi_set_lookup(b->lookup, e, tag);
    ++(m->n[dim]);
  }
  return gmi_from_agm(e);
}

static struct gmi_model_ops ops = {
  .begin   = gmi_base_begin,
  .next    = gmi_base_next,
  .end     = gmi_base_end,
  .dim     = gmi_base_dim,
  .tag     = gmi_base_tag,
  .find    = gmi_null_find,
  .destroy = gmi_base_destroy
};

static struct gmi_model* create(const char* filename)
{
  (void)filename;
  struct gmi_base* m;
  m = calloc(1, sizeof(*m));
  m->model.ops = &ops;
  gmi_base_init(m);
  return &m->model;
}

void gmi_register_null(void)
{
  gmi_register(create, "null");
}
