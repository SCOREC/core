/****************************************************************************** 

  Copyright 2014 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/
#include "gmi_mesh.h"
#include <stdlib.h>

static struct gmi_model_ops ops = {
  .begin    = gmi_base_begin,
  .next     = gmi_base_next,
  .end      = gmi_base_end,
  .dim      = gmi_base_dim,
  .tag      = gmi_base_tag,
  .find     = gmi_base_find,
  .adjacent = gmi_base_adjacent,
  .destroy  = gmi_base_destroy
};

static struct gmi_model* create(const char* filename)
{
  struct gmi_base* m;
  FILE* f;
  f = fopen(filename, "r");
  if (!f)
    gmi_fail("could not open model file");
  m = malloc(sizeof(*m));
  m->model.ops = &ops;
  gmi_base_read_dmg(m, f);
  fclose(f);
  return &m->model;
}

void gmi_register_mesh(void)
{
  gmi_register(create, "dmg");
}
