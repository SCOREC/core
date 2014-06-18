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

static void eval(struct gmi_model* m, struct gmi_ent* e,
      double const p[2], double x[3])
{
  int dim = gmi_dim(m, e);
  if (dim == 2) {
    GF_point((pGFace)e, p, x);
    return;
  }
  if (dim == 1) {
    GE_point((pGEdge)e, p[0], x);
    return;
  }
  if (dim == 0) {
    GV_point((pGVertex)e, x);
  }
  gmi_fail("bad dimension in gmi_sim eval");
}

static void reparam(struct gmi_model* m, struct gmi_ent* from,
      double const from_p[2], struct gmi_ent* to, double to_p[2])
{
  int from_dim, to_dim;
  from_dim = gmi_dim(m, from);
  to_dim = gmi_dim(m, to);
  if ((from_dim == 1) && (to_dim == 2)) {
    GF_edgeReparam((pGFace)to, (pGEdge)from, from_p[0], 1, to_p);
    return;
  }
  if ((from_dim == 0) && (to_dim == 2)) {
    GF_vertexReparam((pGFace)to, (pGVertex)from, to_p);
    return;
  }
  if ((from_dim == 0) && (to_dim == 1)) {
    to_p[0] = GE_vertexReparam((pGEdge)to, (pGVertex)from);
    return;
  }
  gmi_fail("bad dimensions in gmi_sim reparam");
}

static int periodic(struct gmi_model* m, struct gmi_ent* e, int dim)
{
  int md = gmi_dim(m, e);
  if (md == 2)
    return GF_isSurfacePeriodic((pGFace)e, dim);
  if (md == 1)
    return GE_periodic((pGEdge)e);
  return 0;
}

static void range(struct gmi_model* m, struct gmi_ent* e, int dim,
    double r[2])
{
  int md = gmi_dim(m, e);
  if (md == 2)
    return GF_parRange((pGFace)e, dim, &r[0], &r[1]);
  if (md == 1)
    return GE_parRange((pGEdge)e, &r[0], &r[1]);
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
  ops.eval = eval;
  ops.reparam = reparam;
  ops.periodic = periodic;
  ops.range = range;
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
