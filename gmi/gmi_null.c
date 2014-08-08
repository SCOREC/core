/****************************************************************************** 

  Copyright 2014 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/
#include "gmi_null.h"
#include "gmi_mesh.h"
#include <stdlib.h>

struct gmi_ent* gmi_null_find(struct gmi_model* m, int dim, int tag)
{
  int i;
  struct gmi_mesh* m2 = (struct gmi_mesh*)m;
  for (i = 0; i < m->n[dim]; ++i)
    if (m2->tags[dim][i] == tag)
      break;
  if (i == m->n[dim]) {
    ++(m->n[dim]);
    m2->tags[dim] = realloc(m2->tags[dim], m->n[dim] * sizeof(int));
    m2->tags[dim][i] = tag;
  }
  return gmi_identify(dim, i);
}

/* this is different from gmi_mesh_destroy
   because the null model doesn't develop
   a topology, and gmi_mesh_destroy tries
   to free the topology */
static void destroy(struct gmi_model* m)
{
  struct gmi_mesh* mm = (struct gmi_mesh*)m;
  int i;
  for (i = 0; i < 4; ++i)
    free(mm->tags[i]);
  free(mm);
}

static struct gmi_model_ops ops = {
  .begin   = gmi_mesh_begin,
  .next    = gmi_mesh_next,
  .end     = gmi_mesh_end,
  .dim     = gmi_mesh_dim,
  .tag     = gmi_mesh_tag,
  .find    = gmi_null_find,
  .destroy = destroy
};

static struct gmi_model* create(const char* filename)
{
  struct gmi_mesh* m;
  m = calloc(1, sizeof(*m));
  m->model.ops = &ops;
  return &m->model;
}

void gmi_register_null(void)
{
  gmi_register(create, "null");
}
