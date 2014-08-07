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
#include <stdarg.h>
#include <assert.h>

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

int gmi_index(struct gmi_ent* e)
{
  return INDEX(from_ptr(e));
}

struct gmi_ent* gmi_identify(int dim, int idx)
{
  return to_ptr(ENT(dim, idx));
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
  int* i;
  i = malloc(sizeof(int));
  *i = first(to_mesh(m), dim);
  return (struct gmi_iter*)i;
}

struct gmi_ent* gmi_mesh_next(struct gmi_model* m, struct gmi_iter* it)
{
  int* i = (int*)it;
  int e = *i;
  *i = next(to_mesh(m), *i);
  return to_ptr(e);
}

void gmi_mesh_end(struct gmi_model* m, struct gmi_iter* i)
{
  free(i);
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
  index = gmi_index(e);
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
  for (i = 0; i < 4; ++i)
    free(mm->tags[i]);
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
  free(mm);
}

static struct gmi_model_ops mesh_ops = {
  .begin    = gmi_mesh_begin,
  .next     = gmi_mesh_next,
  .end      = gmi_mesh_end,
  .dim      = gmi_mesh_dim,
  .tag      = gmi_mesh_tag,
  .find     = gmi_mesh_find,
  .adjacent = gmi_mesh_adjacent,
  .destroy  = gmi_mesh_destroy
};

void gmi_mesh_create(struct gmi_mesh* m, int n[4])
{
  int i, j;
  m->model.ops = &mesh_ops;
  memcpy(m->model.n, n, sizeof(m->model.n));
  for (i = 0; i < 4; ++i)
    m->tags[i] = malloc(n[i] * sizeof(int));
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

static void sort_dim(struct gmi_mesh* m, int dim)
{
  qsort(m->tags[dim], m->model.n[dim], sizeof(int), comp_ints);
}

static void my_fscanf(FILE* f, int n, const char* format, ...)
{
  va_list ap;
  int r;
  va_start(ap, format);
  r = vfscanf(f, format, ap);
  va_end(ap);
  assert(r == n);
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
  dim = gmi_dim(&m->model, e);
  if (is_up)
    ae = gmi_find(&m->model, dim + 1, atag);
  else
    ae = gmi_find(&m->model, dim - 1, atag);
  if (!ae)
    return;
  index = gmi_index(e);
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

void read_tags(struct gmi_mesh* m, FILE* f)
{
  int n[4];
  int i,j,k;
  int loops, shells;
  int faces, edges;
  /* read entity counts */
  my_fscanf(f, 4, "%d %d %d %d", &n[3], &n[2], &n[1], &n[0]);
  gmi_mesh_create(m, n);
  /* bounding box */
  my_fscanf(f, 0, "%*f %*f %*f");
  my_fscanf(f, 0, "%*f %*f %*f");
  /* vertices */
  for (i = 0; i < n[0]; ++i)
    my_fscanf(f, 1, "%d %*f %*f %*f", &m->tags[0][i]);
  sort_dim(m, 0);
  /* edges */
  for (i = 0; i < n[1]; ++i)
    my_fscanf(f, 1, "%d %*d %*d", &m->tags[1][i]);
  sort_dim(m, 1);
  /* faces */
  for (i = 0; i < n[2]; ++i) {
    my_fscanf(f, 2, "%d %d", &m->tags[2][i], &loops);
    for (j = 0; j < loops; ++j) {
      my_fscanf(f, 1, "%d", &edges);
      for (k = 0; k < edges; ++k)
        /* tag, direction */
        my_fscanf(f, 0, "%*d %*d");
    }
  }
  sort_dim(m, 2);
  /* regions */
  for (i = 0; i < n[3]; ++i) {
    my_fscanf(f, 2, "%d %d", &m->tags[3][i], &shells);
    for (j = 0; j < shells; ++j) {
      my_fscanf(f, 1, "%d", &faces);
      for (k = 0; k < faces; ++k)
        /* tag, direction */
        my_fscanf(f, 0, "%*d %*d");
    }
  }
  sort_dim(m, 3);
}

void read_topology(struct gmi_mesh* m, FILE* f)
{
  int i,j,k;
  int loops, shells;
  int faces, edges;
  int tag;
  struct gmi_ent* e;
  /* read entity counts */
  my_fscanf(f, 0, "%*d %*d %*d %*d");
  /* bounding box */
  my_fscanf(f, 0, "%*f %*f %*f");
  my_fscanf(f, 0, "%*f %*f %*f");
  /* vertices */
  for (i = 0; i < m->model.n[0]; ++i)
    my_fscanf(f, 0, "%*d %*f %*f %*f");
  /* edges */
  for (i = 0; i < m->model.n[1]; ++i) {
    my_fscanf(f, 1, "%d", &tag);
    e = gmi_find(&m->model, 1, tag);
    for (j = 0; j < 2; ++j) {
      my_fscanf(f, 1, "%d", &tag);
      add_adjacent(m, e, 0, tag);
    }
  }
  /* faces */
  for (i = 0; i < m->model.n[2]; ++i) {
    my_fscanf(f, 2, "%d %d", &tag, &loops);
    e = gmi_find(&m->model, 2, tag);
    for (j = 0; j < loops; ++j) {
      my_fscanf(f, 1, "%d", &edges);
      for (k = 0; k < edges; ++k) {
        /* tag, direction */
        my_fscanf(f, 1, "%d %*d", &tag);
        add_adjacent(m, e, 0, tag);
      }
    }
  }
  /* regions */
  for (i = 0; i < m->model.n[3]; ++i) {
    my_fscanf(f, 2, "%d %d", &tag, &shells);
    e = gmi_find(&m->model, 3, tag);
    for (j = 0; j < shells; ++j) {
      my_fscanf(f, 1, "%d", &faces);
      for (k = 0; k < faces; ++k) {
        my_fscanf(f, 1, "%d %*d", &tag);
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
  read_tags(m, f);
  rewind(f);
  read_topology(m, f);
  fclose(f);
}

void gmi_write_dmg(struct gmi_model* m, const char* filename)
{
  struct gmi_iter* it;
  struct gmi_ent* e;
  struct gmi_set* s;
  FILE* f = fopen(filename, "w");
  int i;
  /* entity counts */
  fprintf(f, "%d %d %d %d\n", m->n[3], m->n[2], m->n[1], m->n[0]);
  /* bounding box */
  fprintf(f, "0 0 0\n");
  fprintf(f, "0 0 0\n");
  /* vertices */
  it = gmi_begin(m, 0);
  while ((e = gmi_next(m, it))) {
    fprintf(f, "%d 0 0 0\n", gmi_tag(m, e));
  }
  gmi_end(m, it);
  /* edges */
  it = gmi_begin(m, 1);
  while ((e = gmi_next(m, it))) {
    s = gmi_adjacent(m, e, 0);
    fprintf(f, "%d ", gmi_tag(m, e));
    if (s->n == 0)
      fprintf(f,"-42 -42\n");
    else if (s->n == 1)
      fprintf(f,"%d -42\n", gmi_tag(m, s->e[0]));
    else
      fprintf(f, "%d %d\n",
          gmi_tag(m, s->e[0]), gmi_tag(m, s->e[1]));
    gmi_free_set(s);
  }
  gmi_end(m, it);
  /* faces */
  it = gmi_begin(m, 2);
  while ((e = gmi_next(m, it))) {
    /* we're going to cheat a bit here and
       treat all edges as one giant loop */
    fprintf(f, "%d 1\n", gmi_tag(m, e));
    s = gmi_adjacent(m, e, 1);
    fprintf(f, "%d\n", s->n);
    for (i = 0; i < s->n; ++i)
      fprintf(f, " %d 0\n", gmi_tag(m, s->e[i]));
    gmi_free_set(s);
  }
  gmi_end(m, it);
  /* regions */
  it = gmi_begin(m, 3);
  while ((e = gmi_next(m, it))) {
    /* same sort of cheat, all faces are one shell */
    fprintf(f, "%d 1\n", gmi_tag(m, e));
    s = gmi_adjacent(m, e, 2);
    fprintf(f, "%d\n", s->n);
    for (i = 0; i < s->n; ++i)
      fprintf(f, " %d 0\n", gmi_tag(m, s->e[i]));
    gmi_free_set(s);
  }
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
