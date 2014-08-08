/****************************************************************************** 

  Copyright 2014 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/
#include "gmi_base.h"
#include <stdlib.h>
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

int gmi_base_index(struct gmi_ent* e)
{
  return INDEX(from_ptr(e));
}

struct gmi_ent* gmi_base_identify(int dim, int idx)
{
  return to_ptr(ENT(dim, idx));
}

static struct gmi_base* to_base(struct gmi_model* m)
{
  return (struct gmi_base*)m;
}

static int first(struct gmi_base* m, int dim)
{
  if (m->model.n[dim])
    return ENT(dim, 0);
  return -1;
}

static int next(struct gmi_base* m, int e)
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

struct gmi_iter* gmi_base_begin(struct gmi_model* m, int dim)
{
  int* i;
  i = malloc(sizeof(int));
  *i = first(to_base(m), dim);
  return (struct gmi_iter*)i;
}

struct gmi_ent* gmi_base_next(struct gmi_model* m, struct gmi_iter* it)
{
  int* i = (int*)it;
  int e = *i;
  *i = next(to_base(m), *i);
  return to_ptr(e);
}

void gmi_base_end(struct gmi_model* m, struct gmi_iter* i)
{
  free(i);
}

int gmi_base_dim(struct gmi_model* m, struct gmi_ent* e)
{
  return DIM(from_ptr(e));
}

int gmi_base_tag(struct gmi_model* m, struct gmi_ent* e)
{
  int id = from_ptr(e);
  return to_base(m)->tags[DIM(id)][INDEX(id)];
}

static int comp_ints(const void* p, const void* q)
{
  int const* a = p;
  int const* b = q;
  return *a - *b;
}

struct gmi_ent* gmi_base_find(struct gmi_model* m, int dim, int tag)
{
  struct gmi_base* mm = to_base(m);
  int* found;
  int i;
  found = bsearch(&tag, mm->tags[dim], mm->model.n[dim],
      sizeof(int), comp_ints);
  if (!found)
    return NULL;
  i = found - mm->tags[dim];
  return to_ptr(ENT(dim, i));
}

void gmi_base_destroy(struct gmi_model* m)
{
  struct gmi_base* mm = to_base(m);
  int i;
  for (i = 0; i < 4; ++i)
    free(mm->tags[i]);
  free(mm);
}

static void create(struct gmi_base* m, int n[4])
{
  int i;
  memcpy(m->model.n, n, sizeof(m->model.n));
  for (i = 0; i < 4; ++i)
    m->tags[i] = malloc(n[i] * sizeof(int));
}

static void sort_dim(struct gmi_base* m, int dim)
{
  qsort(m->tags[dim], m->model.n[dim], sizeof(int), comp_ints);
}

void gmi_base_read_dmg(struct gmi_base* m, FILE* f)
{
  int n[4];
  int i,j,k;
  int loops, shells;
  int faces, edges;
  /* read entity counts */
  gmi_fscanf(f, 4, "%d %d %d %d", &n[3], &n[2], &n[1], &n[0]);
  create(m, n);
  /* bounding box */
  gmi_fscanf(f, 0, "%*f %*f %*f");
  gmi_fscanf(f, 0, "%*f %*f %*f");
  /* vertices */
  for (i = 0; i < n[0]; ++i)
    gmi_fscanf(f, 1, "%d %*f %*f %*f", &m->tags[0][i]);
  sort_dim(m, 0);
  /* edges */
  for (i = 0; i < n[1]; ++i)
    gmi_fscanf(f, 1, "%d %*d %*d", &m->tags[1][i]);
  sort_dim(m, 1);
  /* faces */
  for (i = 0; i < n[2]; ++i) {
    gmi_fscanf(f, 2, "%d %d", &m->tags[2][i], &loops);
    for (j = 0; j < loops; ++j) {
      gmi_fscanf(f, 1, "%d", &edges);
      for (k = 0; k < edges; ++k)
        /* tag, direction */
        gmi_fscanf(f, 0, "%*d %*d");
    }
  }
  sort_dim(m, 2);
  /* regions */
  for (i = 0; i < n[3]; ++i) {
    gmi_fscanf(f, 2, "%d %d", &m->tags[3][i], &shells);
    for (j = 0; j < shells; ++j) {
      gmi_fscanf(f, 1, "%d", &faces);
      for (k = 0; k < faces; ++k)
        /* tag, direction */
        gmi_fscanf(f, 0, "%*d %*d");
    }
  }
  sort_dim(m, 3);
}

