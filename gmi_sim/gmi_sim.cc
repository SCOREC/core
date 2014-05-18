/****************************************************************************** 

  Copyright 2014 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/
#include "gmi_sim.h"
#include <gmi.h>
#include <stdlib.h>
#include <gmi.h>
#include <SimModel.h>

struct sim_model {
  struct gmi_model model;
  SGModel* sim;
};

struct sim_iter {
  union {
    GVIter v;
    GEIter e;
    GFIter f;
    GRIter r;
  } i;
  int dim;
};

extern "C" {

static gmi_iter* begin(gmi_model* m, int dim)
{
  sim_model* mm;
  sim_iter* i;
  mm = (sim_model*)m;
  i = (sim_iter*)malloc(sizeof(*i));
  i->dim = dim;
  if (dim == 0)
    i->i.v = GM_vertexIter(mm->sim);
  else if (dim == 1)
    i->i.e = GM_edgeIter(mm->sim);
  else if (dim == 2)
    i->i.f = GM_faceIter(mm->sim);
  else if (dim == 3)
    i->i.r = GM_regionIter(mm->sim);
  return (gmi_iter*)i;
}

static gmi_ent* next(gmi_model* m, gmi_iter* i)
{
  sim_iter* si;
  si = (sim_iter*)i;
  if (si->dim == 0)
    return (gmi_ent*)GVIter_next(si->i.v);
  if (si->dim == 1)
    return (gmi_ent*)GEIter_next(si->i.e);
  if (si->dim == 2)
    return (gmi_ent*)GFIter_next(si->i.f);
  if (si->dim == 3)
    return (gmi_ent*)GRIter_next(si->i.r);
  return 0;
}

static void end(gmi_model* m, gmi_iter* i)
{
  sim_iter* si;
  si = (sim_iter*)i;
  if (si->dim == 0)
    return GVIter_delete(si->i.v);
  if (si->dim == 1)
    return GEIter_delete(si->i.e);
  if (si->dim == 2)
    return GFIter_delete(si->i.f);
  if (si->dim == 3)
    return GRIter_delete(si->i.r);
  free(si);
}

static int get_dim(gmi_model* m, gmi_ent* e)
{
  return GEN_type((pGEntity)e);
}

static int get_tag(gmi_model* m, gmi_ent* e)
{
  return GEN_tag((pGEntity)e);
}

static gmi_ent* find(gmi_model* m, int dim, int tag)
{
  sim_model* mm = (sim_model*)m;
  return (gmi_ent*)GM_entityByTag(mm->sim, dim, tag);
}

static void destroy(gmi_model* m)
{
  sim_model* mm = (sim_model*)m;
  GM_release(mm->sim);
  free(mm);
}

static struct gmi_model_ops ops;

static gmi_model* create(const char* filename)
{
  return gmi_import_sim(GM_load(filename, NULL, NULL));
}

void gmi_register_sim(void)
{
  ops.begin = begin;
  ops.next = next;
  ops.end = end;
  ops.dim = get_dim;
  ops.tag = get_tag;
  ops.find = find;
  ops.destroy = destroy;
  gmi_register(create, "smd");
}

gmi_model* gmi_import_sim(SGModel* sm)
{
  sim_model* m;
  m = (sim_model*)malloc(sizeof(*m));
  m->model.ops = &ops;
  m->sim = sm;
  m->model.n[0] = GM_numVertices(sm);
  m->model.n[1] = GM_numEdges(sm);
  m->model.n[2] = GM_numFaces(sm);
  m->model.n[3] = GM_numRegions(sm);
  return &m->model;
}

}
