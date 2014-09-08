/****************************************************************************** 

  Copyright 2014 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/
#include "gmi_null.h"
#include "gmi_base.h"
#include <stdlib.h>

struct gmi_ent* gmi_null_find(struct gmi_model* m, int dim, int tag)
{
  int i;
  struct gmi_base* m2 = (struct gmi_base*)m;
  for (i = 0; i < m->n[dim]; ++i)
    if (m2->tags[dim][i] == tag)
      break;
  if (i == m->n[dim]) {
    ++(m->n[dim]);
    m2->tags[dim] = realloc(m2->tags[dim], m->n[dim] * sizeof(int));
    m2->tags[dim][i] = tag;
  }
  return gmi_base_identify(dim, i);
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
  return &m->model;
}

void gmi_register_null(void)
{
  gmi_register(create, "null");
}
