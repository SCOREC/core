/****************************************************************************** 

  Copyright 2014 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/

#include "mds_apf.h"
#include <stdlib.h>
#include <assert.h>
#include <PCU.h>

struct mds_apf* mds_apf_create(struct gmi_model* model, int d,
    mds_id cap[MDS_TYPES])
{
  struct mds_apf* m;
  int t;
  m = malloc(sizeof(*m));
  mds_create(&(m->mds),d,cap);
  mds_create_tags(&(m->tags));
  m->point = malloc(cap[MDS_VERTEX] * sizeof(*(m->point)));
  m->param = malloc(cap[MDS_VERTEX] * sizeof(*(m->param)));
  for (t = 0; t < MDS_TYPES; ++t)
    m->model[t] = malloc(cap[t] * sizeof(*(m->model[t])));
  m->user_model = model;
  for (t = 0; t < MDS_TYPES; ++t)
    m->parts[t] = calloc(cap[t], sizeof(*(m->parts[t])));
  mds_create_net(&m->remotes);
  mds_create_net(&m->matches);
  return m;
}

void mds_apf_destroy(struct mds_apf* m)
{
  int t;
  mds_destroy_net(&m->matches, &m->mds);
  mds_destroy_net(&m->remotes, &m->mds);
  for (t = 0; t < MDS_TYPES; ++t)
    free(m->model[t]);
  for (t = 0; t < MDS_TYPES; ++t)
    free(m->parts[t]);
  free(m->point);
  free(m->param);
  mds_destroy_tags(&(m->tags));
  mds_destroy(&(m->mds));
  free(m);
}

double* mds_apf_point(struct mds_apf* m, mds_id e)
{
  return m->point[mds_index(e)];
}

double* mds_apf_param(struct mds_apf* m, mds_id e)
{
  return m->param[mds_index(e)];
}

struct gmi_ent* mds_apf_model(struct mds_apf* m, mds_id e)
{
  return m->model[mds_type(e)][mds_index(e)];
}

mds_id mds_apf_create_entity(
    struct mds_apf* m, int type, struct gmi_ent* model, mds_id* from)
{
  int t;
  mds_id old_cap[MDS_TYPES];
  mds_id e;
  mds_id i;
  old_cap[type] = m->mds.cap[type];
  e = mds_create_entity(&(m->mds),type,from);
  i = mds_index(e);
  if (m->mds.cap[type] != old_cap[type]) {
    for (t = 0; t < MDS_TYPES; ++t)
      if (t != type)
        old_cap[t] = m->mds.cap[t];
    mds_grow_tags(&(m->tags),&(m->mds),old_cap);
    if (type == MDS_VERTEX) {
      m->point = realloc(m->point,m->mds.cap[type] * sizeof(*(m->point)));
      m->param = realloc(m->param,m->mds.cap[type] * sizeof(*(m->param)));
    }
    m->model[type] = realloc(m->model[type],
        m->mds.cap[type] * sizeof(*(m->model[type])));
    m->parts[type] = realloc(m->parts[type],
        m->mds.cap[type] * sizeof(*(m->parts[type])));
    mds_grow_net(&m->remotes, &m->mds, old_cap);
    mds_grow_net(&m->matches, &m->mds, old_cap);
  }
  m->model[type][i] = model;
  m->parts[type][i] = NULL;
  if (type == MDS_VERTEX) {
    m->point[i][0] = m->point[i][1] = m->point[i][2] = 0;
    m->param[i][0] = m->param[i][1] = 0;
  }
  return e;
}

void mds_apf_destroy_entity(struct mds_apf* m, mds_id e)
{
  struct mds_tag* t;
  for (t = m->tags.first; t; t = t->next)
    if (mds_has_tag(t,e))
      mds_take_tag(t,e);
  mds_set_copies(&m->remotes, &m->mds, e, NULL);
  mds_set_copies(&m->matches, &m->mds, e, NULL);
  mds_destroy_entity(&(m->mds),e);
}

void* mds_get_part(struct mds_apf* m, mds_id e)
{
  return m->parts[mds_type(e)][mds_index(e)];
}

void mds_set_part(struct mds_apf* m, mds_id e, void* p)
{
  m->parts[mds_type(e)][mds_index(e)] = p;
}

struct gmi_ent* mds_find_model(struct mds_apf* m, int dim, int id)
{
  return gmi_find(m->user_model, dim, id);
}

int mds_model_dim(struct mds_apf* m, struct gmi_ent* model)
{
  return gmi_dim(m->user_model, model);
}

int mds_model_id(struct mds_apf* m, struct gmi_ent* model)
{
  return gmi_tag(m->user_model, model);
}

static void downs_to_copy(struct mds_set* s,
    struct mds_copy c)
{
  int i;
  PCU_COMM_PACK(c.p, c.e);
  for (i = 0; i < s->n; ++i)
    PCU_COMM_PACK(c.p, s->e[i]);
}

static void downs_to_copies(
    struct mds* m, mds_id e, struct mds_copies* c)
{
  int i;
  struct mds_set s;
  mds_get_adjacent(m, e, mds_dim[mds_type(e)] - 1, &s);
  for (i = 0; i < c->n; ++i)
    downs_to_copy(&s, c->c[i]);
}

static void change_down(struct mds* m, mds_id e, struct mds_set* s)
{
  mds_id e2;
  /* note: this sortof hack relies on the LIFO property
     of create/destroy */
  mds_destroy_entity(m, e);
  e2 = mds_create_entity(m, mds_type(e), s->e);
  assert(e2 == e);
}

static int has_copy(struct mds_net* net, mds_id e, struct mds_copy c)
{
  int i;
  struct mds_copies* cs;
  cs = mds_get_copies(net, e);
  assert(cs);
  for (i = 0; i < cs->n; ++i)
    if ((cs->c[i].p == c.p)&&(cs->c[i].e == c.e))
      return 1;
  return 0;
}

static int compare_copy_sets(struct mds_net* net,
    struct mds_set* local,
    int from,
    struct mds_set* remote)
{
  int i;
  struct mds_copy c;
  c.p = from;
  assert(local->n == remote->n);
  for (i = 0; i < local->n; ++i) {
    c.e = remote->e[i];
    if (!has_copy(net, local->e[i], c))
      return 0;
  }
  return 1;
}

/* accepts negative offsets, which are of opposite curl */
static void rotate_set(struct mds_set* in, int r, struct mds_set* out)
{
  int i;
  out->n = in->n;
  if (r < 0)
    for (i = 0; i < in->n; ++i)
      out->e[i] = in->e[(in->n - r - i) % in->n];
  else
    for (i = 0; i < in->n; ++i)
      out->e[i] = in->e[(i + r) % in->n];
}

static int recv_down_copies(struct mds_net* net, struct mds* m)
{
  mds_id e;
  struct mds_set s;
  struct mds_set rs;
  struct mds_set s2;
  int i;
  int from = PCU_Comm_Sender();
  PCU_COMM_UNPACK(e);
  mds_get_adjacent(m, e, mds_dim[mds_type(e)] - 1, &s);
  rs.n = s.n;
  for (i = 0; i < s.n; ++i)
    PCU_COMM_UNPACK(rs.e[i]);
  if (compare_copy_sets(net, &s, from, &rs))
    return 0;
  for (i = -s.n; i < s.n; ++i) {
    rotate_set(&s, i, &s2);
    if (compare_copy_sets(net, &s2, from, &rs)) {
      change_down(m, e, &s2);
      return 1;
    }
  }
  abort();
}

static int copy_less(struct mds_copy a, struct mds_copy b)
{
  if (a.p != b.p)
    return a.p < b.p;
  return a.e < b.e;
}

static int owns_copies(mds_id e, struct mds_copies* c)
{
  int i;
  struct mds_copy mc;
  mc.e = e;
  mc.p = PCU_Comm_Self();
  for (i = 0; i < c->n; ++i)
    if (copy_less(c->c[i], mc))
      return 0;
  return 1;
}

static int align_copies(struct mds_net* net, struct mds* m)
{
  int d;
  mds_id e;
  struct mds_copies* c;
  int did_change = 0;
  PCU_Comm_Begin();
  for (d = 1; d < m->d; ++d)
    for (e = mds_begin(m, d); e != MDS_NONE; e = mds_next(m, e)) {
      c = mds_get_copies(net, e);
      if (!c)
        continue;
      if (owns_copies(e, c))
        downs_to_copies(m, e, c);
    }
  PCU_Comm_Send();
  while (PCU_Comm_Receive())
    if (recv_down_copies(net, m))
      did_change = 1;
  return PCU_Or(did_change);
}

int mds_align_matches(struct mds_apf* m)
{
  return align_copies(&m->matches, &m->mds);
}

int mds_align_remotes(struct mds_apf* m)
{
  return align_copies(&m->remotes, &m->mds);
}

/* uses the null model to classify a mesh
   that did not have a model in such a way
   that verification will accept it */
void mds_derive_model(struct mds_apf* m)
{
  int d;
  mds_id e;
  int i;
  mds_id de;
  struct mds_set s;
  struct mds_copies* c;
  struct gmi_ent* interior = mds_find_model(m, m->mds.d, 0);
  struct gmi_ent* boundary = mds_find_model(m, m->mds.d - 1, 0);
  /* first classify everything to the interior */
  for (d = 0; d <= m->mds.d; ++d)
    for (e = mds_begin(&m->mds, d);
         e != MDS_NONE;
         e = mds_next(&m->mds, e))
      m->model[mds_type(e)][mds_index(e)] = interior;
  /* then if a face has neither two adjacent elements
     nor a remote copy, classify its closure onto the model boundary */
  for (e = mds_begin(&m->mds, m->mds.d - 1);
       e != MDS_NONE;
       e = mds_next(&m->mds, e)) {
    mds_get_adjacent(&m->mds, e, m->mds.d, &s);
    c = mds_get_copies(&m->remotes, e);
    if (c || s.n == 2)
      continue;
    for (d = 0; d < m->mds.d; ++d) {
      mds_get_adjacent(&m->mds, e, d, &s);
      for (i = 0; i < s.n; ++i) {
        de = s.e[i];
        m->model[mds_type(de)][mds_index(de)] = boundary;
      }
    }
  }
}
