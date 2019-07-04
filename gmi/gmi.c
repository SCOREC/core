/****************************************************************************** 

  Copyright 2014 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/
#include "gmi.h"
#include <lionPrint.h>
#include <stdlib.h>
#include <string.h>
#include <pcu_util.h>

struct creator {
  gmi_creator f;
  char* ext;
};

struct creators {
  int n;
  struct creator c[1];
};

static struct creators* ctors = NULL;

struct gmi_set* gmi_make_set(int n)
{
  int extra;
  struct gmi_set* s;
  extra = n - 1;
  /* undefined behavior sanitizer complains if we cut
     below the size in the struct definition */
  if (extra < 1)
    extra = 1;
  s = malloc(sizeof(*s) + extra * sizeof(struct gmi_ent*));
  s->n = n;
  return s;
}

void gmi_free_set(struct gmi_set* s)
{
  free(s);
}

struct gmi_iter* gmi_begin(struct gmi_model* m, int dim)
{
  return m->ops->begin(m, dim);
}

struct gmi_ent* gmi_next(struct gmi_model* m, struct gmi_iter* i)
{
  return m->ops->next(m, i);
}

void gmi_end(struct gmi_model* m, struct gmi_iter* i)
{
  m->ops->end(m, i);
}

int gmi_dim(struct gmi_model* m, struct gmi_ent* e)
{
  return m->ops->dim(m, e);
}

int gmi_tag(struct gmi_model* m, struct gmi_ent* e)
{
  return m->ops->tag(m, e);
}

struct gmi_ent* gmi_find(struct gmi_model* m, int dim, int tag)
{
  return m->ops->find(m, dim, tag);
}

// one-level adjacency only
struct gmi_set* gmi_adjacent(struct gmi_model* m, struct gmi_ent* e, int dim)
{
  if (m->ops->adjacent)
    return m->ops->adjacent(m, e, dim);
  return gmi_make_set(0);
}

int gmi_can_eval(struct gmi_model* m)
{
  return m->ops->eval != NULL;
}

int gmi_can_get_closest_point(struct gmi_model* m)
{
  return m->ops->closest_point != NULL;
}

int gmi_has_normal(struct gmi_model* m)
{
  return m->ops->normal != NULL;
}

void gmi_eval(struct gmi_model* m, struct gmi_ent* e,
    double const p[2], double x[3])
{
  m->ops->eval(m, e, p, x);
}

void gmi_reparam(struct gmi_model* m, struct gmi_ent* from,
    double const from_p[2], struct gmi_ent* to, double to_p[2])
{
  m->ops->reparam(m, from, from_p, to, to_p);
}

int gmi_periodic(struct gmi_model* m, struct gmi_ent* e, int dim)
{
  return m->ops->periodic(m, e, dim);
}

void gmi_range(struct gmi_model* m, struct gmi_ent* e, int dim,
    double r[2])
{
  m->ops->range(m, e, dim, r);
}

void gmi_closest_point(struct gmi_model* m, struct gmi_ent* e,
    double const from[3], double to[3], double to_p[2])
{
  m->ops->closest_point(m, e, from, to, to_p);
}

void gmi_normal(struct gmi_model* m, struct gmi_ent* e,
    double const p[2], double n[3])
{
  m->ops->normal(m, e, p, n);
}

void gmi_first_derivative(struct gmi_model* m, struct gmi_ent* e,
    double const p[2], double t0[3], double t1[3])
{
  m->ops->first_derivative(m, e, p, t0, t1);
}

int gmi_is_point_in_region(struct gmi_model* m, struct gmi_ent* e,
    double point[3])
{
  return m->ops->is_point_in_region(m, e, point);
}

int gmi_is_in_closure_of(struct gmi_model* m, struct gmi_ent* e,
    struct gmi_ent* et)
{
  return m->ops->is_in_closure_of(m, e, et);
}

void gmi_bbox(struct gmi_model* m, struct gmi_ent* e,
    double bmin[3], double bmax[3])
{
  return m->ops->bbox(m, e, bmin, bmax);
}

int gmi_is_discrete_ent(struct gmi_model* m, struct gmi_ent* e)
{
  return m->ops->is_discrete_ent(m, e);
}

void gmi_destroy(struct gmi_model* m)
{
  m->ops->destroy(m);
}

static void free_ctors(void)
{
  int i;
  for (i = 0; i < ctors->n; ++i)
    free(ctors->c[i].ext);
  free(ctors);
}

static char* copy_string(const char* a)
{
  char* b;
  b = malloc(strlen(a) + 1);
  strcpy(b, a);
  return b;
}

void gmi_register(gmi_creator f, const char* ext)
{
  if (!ctors) {
    ctors = malloc(sizeof(struct creators));
    ctors->n = 1;
    ctors->c[0].f = f;
    ctors->c[0].ext = copy_string(ext);
    atexit(free_ctors);
  } else {
    ++ctors->n;
    ctors = realloc(ctors, sizeof(struct creators) +
        (ctors->n - 1) * sizeof(struct creator));
    ctors->c[ctors->n - 1].f = f;
    ctors->c[ctors->n - 1].ext = copy_string(ext);
  }
}

int gmi_has_ext(const char* filename, const char* ext)
{
  const char* c = strrchr(filename, '.');
  if (!c)
    gmi_fail("model file name with no extension");
  ++c; /* exclude the dot itself */
  if (!strcmp(c, ext))
    return 1;
  return 0;
}

struct gmi_model* gmi_load(const char* filename)
{
  int i;
  if (!ctors)
    gmi_fail("no models registered before gmi_load");
  for (i = 0; i < ctors->n; ++i)
    if (gmi_has_ext(filename, ctors->c[i].ext))
      return ctors->c[i].f(filename);
  gmi_fail("model file extension not registered");
  return NULL;
}

void gmi_fail(const char* why)
{
  lion_eprint(1,"gmi failed: %s\n",why);
  abort();
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
