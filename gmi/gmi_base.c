/****************************************************************************** 

  Copyright 2014 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/
#include "gmi_base.h"
#include "gmi_lookup.h"
#include "agm.h"
#include <stdlib.h>
#include <string.h>

struct gmi_ent* gmi_from_agm(struct agm_ent e)
{
  char* p = 0;
  int uid;
  if (agm_ent_null(e))
    uid = 0;
  else
    uid = e.id * AGM_ENT_TYPES + e.type + 1;
  p += uid;
  return (struct gmi_ent*)p;
}

struct agm_ent agm_from_gmi(struct gmi_ent* e)
{
  char* p;
  int uid;
  struct agm_ent a;
  p = (char*)e;
  uid = p - ((char*)0);
  if (uid == 0) {
    a.type = 0;
    a.id = -1;
  } else {
    uid -= 1;
    a.type = uid % AGM_ENT_TYPES;
    a.id = uid / AGM_ENT_TYPES;
  }
  return a;
}

int gmi_base_index(struct gmi_ent* e)
{
  return agm_from_gmi(e).id;
}

struct gmi_ent* gmi_base_identify(int dim, int idx)
{
  struct agm_ent e;
  e.type = agm_type_from_dim(dim);
  e.id = idx;
  return gmi_from_agm(e);
}

static struct gmi_base* to_base(struct gmi_model* m)
{
  return (struct gmi_base*)m;
}

struct agm* gmi_base_topo(struct gmi_model* m)
{
  return to_base(m)->topo;
}

struct gmi_iter* gmi_base_begin(struct gmi_model* m, int dim)
{
  struct agm_ent* i;
  i = malloc(sizeof(*i));
  *i = agm_first_ent(to_base(m)->topo, agm_type_from_dim(dim));
  return (struct gmi_iter*)i;
}

struct gmi_ent* gmi_base_next(struct gmi_model* m, struct gmi_iter* it)
{
  struct agm_ent* i;
  struct agm_ent e;
  i = (struct agm_ent*)it;
  e = *i;
  *i = agm_next_ent(to_base(m)->topo, e);
  return gmi_from_agm(e);
}

void gmi_base_end(struct gmi_model* m, struct gmi_iter* i)
{
  (void)m;
  free(i);
}

int gmi_base_dim(struct gmi_model* m, struct gmi_ent* e)
{
  (void)m;
  return agm_dim_from_type(agm_from_gmi(e).type);
}

int gmi_base_tag(struct gmi_model* m, struct gmi_ent* e)
{
  return gmi_get_lookup(to_base(m)->lookup, agm_from_gmi(e));
}

struct gmi_ent* gmi_base_find(struct gmi_model* m, int dim, int tag)
{
  return gmi_from_agm(gmi_look_up(to_base(m)->lookup, agm_type_from_dim(dim), tag));
}

static struct gmi_set* get_up(struct agm* topo, struct agm_ent e)
{
  int n;
  struct gmi_set* s;
  int i;
  struct agm_use u;
  n = agm_use_count_of(topo, e);
  s = gmi_make_set(n);
  i = 0;
  for (u = agm_first_use_of(topo, e);
       !agm_use_null(u);
       u = agm_next_use_of(topo, u)) {
    s->e[i] = gmi_from_agm(agm_bounds(topo, agm_user(topo, u)));
    ++i;
  }
  return s;
}

static struct gmi_set* get_down(struct agm* topo, struct agm_ent e)
{
  int n;
  struct gmi_set* s;
  int i;
  struct agm_bdry b;
  struct agm_use u;
  n = agm_down_count(topo, e);
  s = gmi_make_set(n);
  i = 0;
  for (b = agm_first_bdry_of(topo, e);
       !agm_bdry_null(b);
       b = agm_next_bdry_of(topo, b)) {
    for (u = agm_first_use_by(topo, b);
         !agm_use_null(u);
         u = agm_next_use_by(topo, u)) {
      s->e[i] = gmi_from_agm(agm_used(topo, u));
      ++i;
    }
  }
  return s;
}

struct gmi_set* gmi_base_adjacent(struct gmi_model* m, struct gmi_ent* e,
    int dim)
{
  int from_dim;
  struct agm_ent a;
  a = agm_from_gmi(e);
  from_dim = agm_dim_from_type(a.type);
  if (dim == from_dim - 1)
    return get_down(to_base(m)->topo, a);
  else if (dim == from_dim + 1)
    return get_up(to_base(m)->topo, a);
  gmi_fail("only one-level adjacencies supported");
  return 0;
}

void gmi_base_destroy(struct gmi_model* m)
{
  struct gmi_base* b;
  b = to_base(m);
  gmi_free_lookup(b->lookup);
  agm_free(b->topo);
  free(b);
}

void gmi_base_init(struct gmi_base* m)
{
  m->topo = agm_new();
  m->lookup = gmi_new_lookup(m->topo);
}

void gmi_base_reserve(struct gmi_base* m, int dim, int n)
{
  agm_reserve(m->topo, agm_type_from_dim(dim), n);
  m->model.n[dim] = n;
}

void gmi_base_freeze(struct gmi_model* m)
{
  struct gmi_base* b;
  int i;
  b = to_base(m);
  for (i = 0; i <= 3; ++i)
    gmi_freeze_lookup(b->lookup, i);
}

void gmi_base_unfreeze(struct gmi_model* m)
{
  struct gmi_base* b;
  b = to_base(m);
  gmi_unfreeze_lookups(b->lookup);
}

void gmi_base_set_tag(struct gmi_model* m, struct gmi_ent* e, int tag)
{
  struct gmi_base* b;
  b = to_base(m);
  gmi_set_lookup(b->lookup, agm_from_gmi(e), tag);
}
