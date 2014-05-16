#include "gmi_mesh.h"
#include <stdlib.h>
#include <stddef.h>

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
