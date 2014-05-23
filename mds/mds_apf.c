/****************************************************************************** 

  Copyright 2014 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/

#include "mds_apf.h"
#include <stdlib.h>

struct mds_apf* mds_apf_create(struct gmi_model* model, int d, int cap[MDS_TYPES])
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
  old_cap[type] = m->mds.cap[type];
  e = mds_create_entity(&(m->mds),type,from);
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
  m->model[type][mds_index(e)] = model;
  m->parts[type][mds_index(e)] = NULL;
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
