/****************************************************************************** 

  Copyright 2014 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/
#include "gmi_mesh.h"
#include <stdlib.h>

struct gmi_mesh {
  struct gmi_base base;
  struct gmi_set** up[4];
  struct gmi_set** down[4];
};

static struct gmi_mesh* to_mesh(struct gmi_model* m)
{
  return (struct gmi_mesh*)m;
}

static struct gmi_set* copy_set(struct gmi_set* s)
{
  int i;
  struct gmi_set* s2 = gmi_make_set(s->n);
  for (i = 0; i < s->n; ++i)
    s2->e[i] = s->e[i];
  return s2;
}

struct gmi_set* gmi_mesh_adjacent(struct gmi_model* m,
    struct gmi_ent* e, int dim)
{
  int index, edim;
  struct gmi_mesh* mm = to_mesh(m);
  index = gmi_base_index(e);
  edim = gmi_dim(m, e);
  if (dim == edim - 1) {
    return copy_set(mm->down[edim][index]);
  } else if (dim == edim + 1)
    return copy_set(mm->up[edim][index]);
  else
    gmi_fail("only one-level adjacencies available");
  return NULL;
}

void gmi_mesh_destroy(struct gmi_model* m)
{
  struct gmi_mesh* mm = to_mesh(m);
  int i, j;
  for (i = 1; i < 4; ++i) {
    for (j = 0; j < m->n[i]; ++j)
      gmi_free_set(mm->down[i][j]);
    free(mm->down[i]);
  }
  for (i = 0; i < 3; ++i) {
    for (j = 0; j < m->n[i]; ++j)
      gmi_free_set(mm->up[i][j]);
    free(mm->up[i]);
  }
  gmi_base_destroy(m);
}

static struct gmi_model_ops ops = {
  .begin    = gmi_base_begin,
  .next     = gmi_base_next,
  .end      = gmi_base_end,
  .dim      = gmi_base_dim,
  .tag      = gmi_base_tag,
  .find     = gmi_base_find,
  .adjacent = gmi_mesh_adjacent,
  .destroy  = gmi_mesh_destroy
};

void alloc_topology(struct gmi_mesh* m)
{
  int i, j;
  int* n;
  m->base.model.ops = &ops;
  n = m->base.model.n;
  for (i = 1; i < 4; ++i) {
    m->down[i] = calloc(n[i], sizeof(struct gmi_set*));
    for (j = 0; j < n[i]; ++j)
      m->down[i][j] = gmi_make_set(0);
  }
  m->down[0] = NULL;
  for (i = 0; i < 3; ++i) {
    m->up[i] = calloc(n[i], sizeof(struct gmi_set*));
    for (j = 0; j < n[i]; ++j)
      m->up[i][j] = gmi_make_set(0);
  }
  m->up[3] = NULL;
}

static void add_adjacent(struct gmi_mesh* m, struct gmi_ent* e,
    int is_up, int atag)
{
  int dim;
  struct gmi_ent* ae;
  int index;
  struct gmi_set** ps;
  struct gmi_set* s;
  int i;
  dim = gmi_dim(&m->base.model, e);
  if (is_up)
    ae = gmi_find(&m->base.model, dim + 1, atag);
  else
    ae = gmi_find(&m->base.model, dim - 1, atag);
  if (!ae)
    return;
  index = gmi_base_index(e);
  if (is_up)
    ps = m->up[dim] + index;
  else
    ps = m->down[dim] + index;
  s = gmi_make_set((*ps)->n + 1);
  for (i = 0; i < (*ps)->n; ++i)
    s->e[i] = (*ps)->e[i];
  s->e[(*ps)->n] = ae;
  gmi_free_set(*ps);
  *ps = s;
}

void read_topology(struct gmi_mesh* m, FILE* f)
{
  int i,j,k;
  int loops, shells;
  int faces, edges;
  int tag;
  struct gmi_ent* e;
  struct gmi_model* mod = &m->base.model;
  /* read entity counts */
  gmi_fscanf(f, 0, "%*d %*d %*d %*d");
  /* bounding box */
  gmi_fscanf(f, 0, "%*f %*f %*f");
  gmi_fscanf(f, 0, "%*f %*f %*f");
  /* vertices */
  for (i = 0; i < mod->n[0]; ++i)
    gmi_fscanf(f, 0, "%*d %*f %*f %*f");
  /* edges */
  for (i = 0; i < mod->n[1]; ++i) {
    gmi_fscanf(f, 1, "%d", &tag);
    e = gmi_find(mod, 1, tag);
    for (j = 0; j < 2; ++j) {
      gmi_fscanf(f, 1, "%d", &tag);
      add_adjacent(m, e, 0, tag);
    }
  }
  /* faces */
  for (i = 0; i < mod->n[2]; ++i) {
    gmi_fscanf(f, 2, "%d %d", &tag, &loops);
    e = gmi_find(mod, 2, tag);
    for (j = 0; j < loops; ++j) {
      gmi_fscanf(f, 1, "%d", &edges);
      for (k = 0; k < edges; ++k) {
        /* tag, direction */
        gmi_fscanf(f, 1, "%d %*d", &tag);
        add_adjacent(m, e, 0, tag);
      }
    }
  }
  /* regions */
  for (i = 0; i < mod->n[3]; ++i) {
    gmi_fscanf(f, 2, "%d %d", &tag, &shells);
    e = gmi_find(mod, 3, tag);
    for (j = 0; j < shells; ++j) {
      gmi_fscanf(f, 1, "%d", &faces);
      for (k = 0; k < faces; ++k) {
        gmi_fscanf(f, 1, "%d %*d", &tag);
        add_adjacent(m, e, 0, tag);
      }
    }
  }
}

void gmi_read_dmg(struct gmi_mesh* m, const char* filename)
{
  FILE* f = fopen(filename, "r");
  if (!f)
    gmi_fail("could not open model file");
  gmi_base_read_dmg(&m->base, f);
  alloc_topology(m);
  rewind(f);
  read_topology(m, f);
  fclose(f);
}

static struct gmi_model* create(const char* filename)
{
  struct gmi_mesh* m;
  m = malloc(sizeof(*m));
  gmi_read_dmg(m, filename);
  return &m->base.model;
}

void gmi_register_mesh(void)
{
  gmi_register(create, "dmg");
}
