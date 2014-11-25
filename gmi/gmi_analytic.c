/*****************************************************************************

  Copyright 2014 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

******************************************************************************/
#include "gmi_analytic.h"
#include "gmi_null.h"
#include "gmi_base.h"
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

typedef uint8_t periodic_t[2];
typedef double ranges_t[2][2];

struct gmi_analytic
{
  struct gmi_base base;
  struct agm_tag* f;
  struct agm_tag* periodic;
  struct agm_tag* ranges;
  struct agm_tag* data;
  struct agm_tag* reparam;
  struct agm_tag* reparam_data;
};

static gmi_analytic_fun* f_of(struct gmi_analytic* m, struct agm_ent e)
{
  return agm_tag_at(m->f, AGM_ENTITY, e.type, e.id);
}

static periodic_t* periodic_of(struct gmi_analytic* m, struct agm_ent e)
{
  return agm_tag_at(m->periodic, AGM_ENTITY, e.type, e.id);
}

static ranges_t* ranges_of(struct gmi_analytic* m, struct agm_ent e)
{
  return agm_tag_at(m->ranges, AGM_ENTITY, e.type, e.id);
}

static void** data_of(struct gmi_analytic* m, struct agm_ent e)
{
  return agm_tag_at(m->data, AGM_ENTITY, e.type, e.id);
}

static gmi_reparam_fun* reparam_of(struct gmi_analytic* m, struct agm_use u)
{
  return agm_tag_at(m->reparam, AGM_USE, u.type, u.id);
}

static void** reparam_data_of(struct gmi_analytic* m, struct agm_use u)
{
  return agm_tag_at(m->reparam_data, AGM_USE, u.type, u.id);
}

static struct gmi_analytic* to_model(struct gmi_model* m)
{
  return (struct gmi_analytic*)m;
}

struct gmi_ent* gmi_add_analytic(struct gmi_model* m, int dim, int tag,
    gmi_analytic_fun f, int* periodic, double (*ranges)[2], void* user_data)
{
  struct gmi_analytic* m2;
  struct agm_ent e;
  periodic_t* pp;
  ranges_t* rp;
  int i;
  m2 = to_model(m);
  e = agm_from_gmi(gmi_null_find(m, dim, tag));
  *(f_of(m2, e)) = f;
  pp = periodic_of(m2, e);
  rp = ranges_of(m2, e);
  for (i = 0; i < dim; ++i) {
    (*pp)[i] = periodic[i];
    (*rp)[i][0] = ranges[i][0];
    (*rp)[i][1] = ranges[i][1];
  }
  for (; i < 2; ++i) {
    (*pp)[i] = 0;
    (*rp)[i][0] = 0;
    (*rp)[i][1] = 0;
  }
  *(data_of(m2, e)) = user_data;
  return gmi_from_agm(e);
}

static void eval(struct gmi_model* m, struct gmi_ent* e,
      double const p[2], double x[3])
{
  struct gmi_analytic* m2;
  struct agm_ent a;
  void* u;
  gmi_analytic_fun f;
  m2 = to_model(m);
  a = agm_from_gmi(e);
  u = *(data_of(m2, a));
  f = *(f_of(m2, a));
  (*f)(p, x, u);
}

static void reparam_across(struct gmi_analytic* m, struct agm_use u,
    double const from_p[2], double to_p[2])
{
  void* d;
  gmi_reparam_fun f;
  d = *(reparam_data_of(m, u));
  f = *(reparam_of(m, u));
  (*f)(from_p, to_p, d);
}

static void reparam_path(struct gmi_analytic* m, struct agm_use* path,
    int pathlen, double const from_p[2], double to_p[2])
{
  double tmp[2];
  if (!pathlen) {
    to_p[0] = from_p[0];
    to_p[1] = from_p[1];
    return;
  }
  reparam_across(m, *path, from_p, tmp);
  reparam_path(m, path + 1, pathlen - 1, tmp, to_p);
}

static void reparam(struct gmi_model* m, struct gmi_ent* from,
      double const from_p[2], struct gmi_ent* to, double to_p[2])
{
  struct gmi_analytic* m2;
  struct agm_ent a;
  struct agm_ent b;
  struct agm_use path[4];
  int pathlen;
  m2 = to_model(m);
  a = agm_from_gmi(from);
  b = agm_from_gmi(to);
  pathlen = agm_find_path(m2->base.topo, a, b, path);
  if (pathlen == -1)
    gmi_fail("analytic reparam can't find topology path");
  reparam_path(m2, path, pathlen, from_p, to_p);
}

static int periodic(struct gmi_model* m, struct gmi_ent* e, int dim)
{
  struct gmi_analytic* m2 = to_model(m);
  periodic_t* pp = periodic_of(m2, agm_from_gmi(e));
  return (*pp)[dim];
}

static void range(struct gmi_model* m, struct gmi_ent* e, int dim, double r[2])
{
  struct gmi_analytic* m2 = to_model(m);
  ranges_t* rp = ranges_of(m2, agm_from_gmi(e));
  r[0] = (*rp)[dim][0];
  r[1] = (*rp)[dim][1];
}

static struct gmi_model_ops ops = {
  .begin    = gmi_base_begin,
  .next     = gmi_base_next,
  .end      = gmi_base_end,
  .dim      = gmi_base_dim,
  .tag      = gmi_base_tag,
  .find     = gmi_base_find,
  .adjacent = gmi_base_adjacent,
  .eval     = eval,
  .reparam  = reparam,
  .periodic = periodic,
  .range    = range,
  .destroy  = gmi_base_destroy
};

struct gmi_model* gmi_make_analytic(void)
{
  struct gmi_analytic* m;
  m = calloc(1, sizeof(*m));
  m->base.model.ops = &ops;
  gmi_base_init(&m->base);
  m->f = agm_new_tag(m->base.topo, sizeof(gmi_analytic_fun));
  m->periodic = agm_new_tag(m->base.topo, sizeof(periodic_t));
  m->ranges = agm_new_tag(m->base.topo, sizeof(ranges_t));
  m->data = agm_new_tag(m->base.topo, sizeof(void*));
  m->reparam = agm_new_tag(m->base.topo, sizeof(gmi_reparam_fun));
  m->reparam_data = agm_new_tag(m->base.topo, sizeof(void*));
  return &m->base.model;
}

void* gmi_analytic_data(struct gmi_model* m, struct gmi_ent* e)
{
  struct gmi_analytic* m2 = to_model(m);
  return *(data_of(m2, agm_from_gmi(e)));
}

void gmi_add_analytic_reparam(struct gmi_model* m, struct agm_use u,
    gmi_reparam_fun f, void* user_data)
{
  struct gmi_analytic* m2 = to_model(m);
  *(reparam_of(m2, u)) = f;
  *(reparam_data_of(m2, u)) = user_data;
}

void* gmi_analytic_reparam_data(struct gmi_model* m, struct agm_use u)
{
  struct gmi_analytic* m2 = to_model(m);
  return *(reparam_data_of(m2, u));
}

void gmi_add_analytic_cell(struct gmi_model* m, int dim, int tag)
{
  gmi_null_find(m, dim, tag);
}

void gmi_add_analytic_region(struct gmi_model* m, int tag)
{
  gmi_add_analytic_cell(m, 3, tag);
}
