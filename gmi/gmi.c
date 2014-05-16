/****************************************************************************** 

  (c) 2011-2014 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/
#include "gmi.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

struct creator {
  gmi_creator f;
  char* ext;
};

struct creators {
  int n;
  struct creator c[];
};

static struct creators* ctors = NULL;

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

int gmi_can_eval(struct gmi_model* m)
{
  return m->ops->eval != NULL;
}

void gmi_eval(struct gmi_model* m, struct gmi_ent* e,
    double const p[2], double x[3])
{
  m->ops->eval(m, e, p, x);
}

void gmi_eval_du(struct gmi_model* m, struct gmi_ent* e,
    double const p[2], double t[3])
{
  m->ops->eval_du(m, e, p, t);
}

void gmi_eval_dv(struct gmi_model* m, struct gmi_ent* e,
    double const p[2], double t[3])
{
  m->ops->eval_dv(m, e, p, t);
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
  return m->ops->range(m, e, dim, r);
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
    ctors = malloc(sizeof(struct creators) + sizeof(struct creator));
    ctors->n = 1;
    ctors->c[0].f = f;
    ctors->c[0].ext = copy_string(ext);
    atexit(free_ctors);
  } else {
    ++ctors->n;
    ctors = realloc(ctors, sizeof(struct creators) +
        ctors->n * sizeof(struct creator));
    ctors->c[ctors->n - 1].f = f;
    ctors->c[ctors->n - 1].ext = copy_string(ext);
  }
}

struct gmi_model* gmi_load(const char* filename)
{
  int i;
  const char* ext = strrchr(filename, '.');
  if (!ext)
    gmi_fail("model file name with no extension");
  ++ext; /* exclude the dot itself */
  if (!ctors)
    gmi_fail("no models registered before gmi_load");
  for (i = 0; i < ctors->n; ++i) {
    if (!strcmp(ext, ctors->c[i].ext))
      return ctors->c[i].f(filename);
  }
  gmi_fail("model file extension not registered");
  return NULL;
}

void gmi_fail(const char* why)
{
  fprintf(stderr,"gmi failed: %s\n",why);
  abort();
}
