/****************************************************************************** 

  Copyright 2014 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/
#include "gmi_mesh.h"
#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <string.h>

#define DIM(e) ((e) % 4)
#define INDEX(e) ((e) / 4)
#define ENT(d,i) ((i) * 4 + (d))

static void* to_ptr(int e)
{
  return ((char*)0) + (e + 1);
}

static int from_ptr(void* e)
{
  return (((char*)e) - ((char*)0)) - 1;
}

static struct gmi_mesh* to_mesh(struct gmi_model* m)
{
  return (struct gmi_mesh*)m;
}

static int first(struct gmi_mesh* m, int dim)
{
  if (m->model.n[dim])
    return ENT(dim, 0);
  return -1;
}

static int next(struct gmi_mesh* m, int e)
{
  int dim;
  int i;
  dim = DIM(e);
  i = INDEX(e);
  ++i;
  if (i == m->model.n[dim])
    return -1;
  return ENT(dim, i);
}

struct gmi_iter* gmi_mesh_begin(struct gmi_model* m, int dim)
{
  return to_ptr(first(to_mesh(m), dim));
}

struct gmi_ent* gmi_mesh_next(struct gmi_model* m, struct gmi_iter* i)
{
  struct gmi_ent* e;
  e = (struct gmi_ent*)i;
  i = to_ptr(next(to_mesh(m), from_ptr(i)));
  return e;
}

void gmi_mesh_end(struct gmi_model* m, struct gmi_iter* i)
{
}

int gmi_mesh_dim(struct gmi_model* m, struct gmi_ent* e)
{
  return DIM(from_ptr(e));
}

int gmi_mesh_tag(struct gmi_model* m, struct gmi_ent* e)
{
  int id = from_ptr(e);
  return to_mesh(m)->tags[DIM(id)][INDEX(id)];
}

static int comp_ints(const void* p, const void* q)
{
  int const* a = p;
  int const* b = q;
  return *a - *b;
}

struct gmi_ent* gmi_mesh_find(struct gmi_model* m, int dim, int tag)
{
  struct gmi_mesh* mm = to_mesh(m);
  int* found;
  int i;
  found = bsearch(&tag, mm->tags[dim], mm->model.n[dim],
      sizeof(int), comp_ints);
  if (!found)
    return NULL;
  i = found - mm->tags[dim];
  return to_ptr(ENT(dim, i));
}

void gmi_mesh_destroy(struct gmi_model* m)
{
  struct gmi_mesh* mm = to_mesh(m);
  int i;
  for (i = 0; i < 4; ++i)
    free(mm->tags[i]);
  free(mm);
}

static struct gmi_model_ops mesh_ops = {
  .begin   = gmi_mesh_begin,
  .next    = gmi_mesh_next,
  .end     = gmi_mesh_end,
  .dim     = gmi_mesh_dim,
  .tag     = gmi_mesh_tag,
  .find    = gmi_mesh_find,
  .destroy = gmi_mesh_destroy
};

void gmi_mesh_create(struct gmi_mesh* m, int n[4])
{
  int i;
  m->model.ops = &mesh_ops;
  memcpy(m->model.n, n, sizeof(m->model.n));
  for (i = 0; i < 4; ++i)
    m->tags[i] = malloc(n[i] * sizeof(int));
}

static void sort_dim(struct gmi_mesh* m, int dim)
{
  qsort(m->tags[dim], m->model.n[dim], sizeof(int), comp_ints);
}

void gmi_read_dmg(struct gmi_mesh* m, const char* filename)
{
  int n[4];
  int i,j,k;
  int loops, shells;
  int faces, edges;
  FILE* f = fopen(filename, "r");
  /* read entity counts */
  fscanf(f, "%d %d %d %d", &n[3], &n[2], &n[1], &n[0]);
  gmi_mesh_create(m, n);
  /* bounding box */
  fscanf(f, "%*f %*f %*f");
  fscanf(f, "%*f %*f %*f");
  /* vertices */
  for (i = 0; i < n[0]; ++i)
    fscanf(f, "%d %*f %*f %*f", &m->tags[0][i]);
  sort_dim(m, 0);
  /* edges */
  for (i = 0; i < n[1]; ++i)
    fscanf(f, "%d %*d %*d", &m->tags[1][i]);
  sort_dim(m, 1);
  /* faces */
  for (i = 0; i < n[2]; ++i) {
    fscanf(f, "%d %d", &m->tags[2][i], &loops);
    for (j = 0; j < loops; ++j) {
      fscanf(f, "%d", &edges);
      for (k = 0; k < edges; ++k)
        /* tag, direction */
        fscanf(f, "%*d %*d");
    }
  }
  sort_dim(m, 2);
  /* regions */
  for (i = 0; i < n[3]; ++i) {
    fscanf(f, "%d %d", &m->tags[3][i], &shells);
    for (j = 0; j < shells; ++j) {
      fscanf(f, "%d", &faces);
      for (k = 0; k < faces; ++k)
        /* tag, direction */
        fscanf(f, "%*d %*d");
    }
  }
  sort_dim(m, 3);
  fclose(f);
}

void gmi_write_dmg(struct gmi_model* m, const char* filename)
{
  struct gmi_iter* it;
  struct gmi_ent* e;
  FILE* f = fopen(filename, "w");
  /* entity counts */
  fprintf(f, "%d %d %d %d\n", m->n[3], m->n[2], m->n[1], m->n[0]);
  /* bounding box */
  fprintf(f, "0 0 0\n");
  fprintf(f, "0 0 0\n");
  /* vertices */
  it = gmi_begin(m, 0);
  while ((e = gmi_next(m, it)))
    fprintf(f, "%d 0 0 0\n", gmi_tag(m, e));
  gmi_end(m, it);
  /* edges */
  it = gmi_begin(m, 1);
  while ((e = gmi_next(m, it)))
    fprintf(f, "%d 0 0\n", gmi_tag(m, e));
  gmi_end(m, it);
  /* faces */
  it = gmi_begin(m, 2);
  while ((e = gmi_next(m, it)))
    fprintf(f, "%d 0\n", gmi_tag(m, e));
  gmi_end(m, it);
  /* regions */
  it = gmi_begin(m, 3);
  while ((e = gmi_next(m, it)))
    fprintf(f, "%d 0\n", gmi_tag(m, e));
  gmi_end(m, it);
  fclose(f);
}

static struct gmi_model* mesh_creator(const char* filename)
{
  struct gmi_mesh* m;
  m = malloc(sizeof(*m));
  gmi_read_dmg(m, filename);
  return &m->model;
}

void gmi_register_mesh(void)
{
  gmi_register(mesh_creator, "dmg");
}
