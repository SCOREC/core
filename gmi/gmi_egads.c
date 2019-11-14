/******************************************************************************

  Copyright 2014 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/
#include "gmi_egads.h"

// initialize to NULL, will be properly set by `gmi_egads_start`
ego *eg_context = NULL;
// initialize to NULL, will be properly set by `gmi_egads_load`
ego *eg_model = NULL;
// initialize to NULL, will be properly set by `gmi_egads_load`
ego *eg_body = NULL;

static struct gmi_iter* begin(struct gmi_model* m, int dim)
{
  ego *eg_ents;
  if (dim = 0)
    EG_getBodyTopos(eg_body, NULL, NODE, NULL, eg_ents);
  else if (dim = 1)
    EG_getBodyTopos(eg_body, NULL, EDGE, NULL, eg_ents);
  else if (dim = 2)
    EG_getBodyTopos(eg_body, NULL, FACE, NULL, eg_ents);
  else if (dim = 3)
    EG_getBodyTopos(eg_body, NULL, SHELL, NULL, eg_ents); // BODY?
  return (struct gmi_iter*)eg_ents;
}

static struct gmi_ent* next(struct gmi_model* m, struct gmi_iter* i)
{
  ego *eg_ents = (ego*)i;
  return ++eg_ents;
}

static void end(struct gmi_model* m, struct gmi_iter* i)
{
  // I think this will create a memory leak as it won't free any of the
  // values that came before
  ego *eg_ents = (ego*)i;
  free(i);
}

static int gmi_dim(struct gmi_model* m, struct gmi_ent* e)
{
  ego *eg_ent = (ego*)e;
  int *ent_type;
  EG_getInfo(eg_ent, ent_type, NULL, NULL, NULL, NULL);
  if (ent_type == NODE)
    return 0;
  else if (ent_type == EDGE)
    return 1;
  else if (ent_type == FACE)
    return 2;
  else if (ent_type == SHELL) // BODY?
    return 3;
  return -1;
}

static int gmi_tag(struct gmi_model* m, struct gmi_ent* e)
{
  /// implement
  // ego *eg_ent = (ego*)e;
}

static struct gmi_ent* find(struct gmi_model* m, int dim, int tag)
{
  /// implement
}

static struct gmi_set* adjacent(struct gmi_model* m, 
                                struct gmi_ent* e, 
                                int dim)
{
  ego *eg_ent = (ego*)e;
  // EG_getBodyTopos(eg_body, e, ...); something like this
}

static void eval(struct gmi_model* m, 
                 struct gmi_ent* e,
                 double const p[2],
                 double x[3])
{
  int *results;
  ego *eg_ent = (ego*)e;
  EG_evaluate(eg_ent, p, results);
  x[0] = results[0];
  x[1] = results[1];
  x[2] = results[2];
  free(results)
}

static void reparam(struct gmi_model* m,
                    struct gmi_ent* from,
                    double const from_p[2],
                    struct gmi_ent* to,
                    double to_p[2])
{
  /// implement
}

static int periodic(struct gmi_model* m,
                    struct gmi_ent* e,
                    int dim)
{
  int ent_dim = gmi_dim(m, e);
  int *periodic;
  ego *eg_ent = (ego*)e;
  EG_getRange(eg_ent, NULL, periodic);

  if (dim == 2) // v direction
  {
    if (ent_dim == 2) // FACE
    {
      if (periodic == 0)
      {
        return 0;
      }
      else if (periodic == 2)
      {
        return 1;
      }
    }
    else
      gmi_fail("v direction only exists for faces");
  }
  return periodic;
}

static void range(struct gmi_model* m,
                  struct gmi_ent* e,
                  int dim,
                  double r[2])
{
  int ent_dim = gmi_dim(m, e);
  int *range;
  ego *eg_ent = (ego*)e;
  EG_getRange(eg_ent, range, NULL);
  if (dim == 2)
  {
    if (ent_dim == 2)
    {
      r[0] = range[2];
      r[1] = range[3];
    }
    else 
      gmi_fail("v direction only exists for faces");
  }
  else if (dim == 1)
  {
    r[0] = range[0];
    r[1] = range[1];
  }
}

static void closest_point(struct gmi_model* m,
                          struct gmi_ent* e, 
                          double const from[3],
                          double to[3],
                          double to_p[2])
{
  ego *eg_ent = (ego*)e;
  EG_invEvaluate(e, &from[0], &to_p[0], &to[0]);
}

static void normal(struct gmi_model* m,
                   struct gmi_ent* e,
                   double const p[2],
                   double n[3])
{
  /// implement
}

static void first_derivative(struct gmi_model* m,
                             struct gmi_ent* e,
                             double const p[2],
                             double t0[3],
                             double t1[3])
{
  int ent_dim = gmi_dim(m, e);
  double *results;
  ego *eg_ent = (ego*)e;
  EG_evaluate(eg_ent, p, results);
  t0[0] = results[3];
  t0[1] = results[4];
  t0[2] = results[5];
  if (ent_dim == 2)
  {
    t1[0] = results[6];
    t1[2] = results[7];
    t1[3] = results[8];
  }
}

static int is_point_in_region(struct gmi_model* m,
                              struct gmi_ent* e,
                              double p[3])
{
  ego *eg_ent = (ego*)e;
  int status = EG_inTopology(e, p);
  if (status == EGADS_SUCCESS)
    return 1;
  else
    return 0;
}

static void bbox(struct gmi_model* m,
                 struct gmi_ent* e,
                 double bmin[3],
                 double bmax[3])
{
  double box[6];
  ego *eg_ent = (ego*)e;
  EG_getBoundingBox(eg_ent, box);
  bmin[0] = box[0];
  bmin[1] = box[1];
  bmin[2] = box[2];
  bmax[0] = box[3];
  bmax[1] = box[4];
  bmax[2] = box[5];
}

static int is_in_closure_of(struct gmi_model* m,
                            struct gmi_ent* e,
                            struct gmi_ent* et)
{
  /// implement
}

static int is_discrete_ent(struct gmi_model* m, struct gmi_ent* e)
{
  /// implement
}

static void destroy(struct gmi_model* m)
{
  /// implement
}

static struct gmi_model_ops ops;

/// TODO: Come up with a better flag? 
// #ifdef HAVE_EGADS
#if 1
static struct gmi_model* gmi_egads_load(const char* filename)
{
  int status = EG_loadModel(*eg_context, 0, filename, eg_model);
  if (status != EGADS_SUCCESS)
  {
    char str[50]; // big enough
    sprintf(str, "EGADS failed to load model with error code: %d", status);
    gmi_fail(str);
  }

  /// TODO: only store the outputs I need, replace the rest with NULL
  int oclass, mtype, nbody, *senses;
  ego geom, *eg_bodies,
  status = EG_getTopology(eg_model, &geom, &oclass, &mtype, NULL, &nbody,
                          &eg_bodies, &senses);
  if (status != EGADS_SUCCESS)
  {
    char str[50]; // big enough
    sprintf(str, "EGADS failed to get bodies with error code: %d", status);
    gmi_fail(str);
  }
  else if (nbody > 1)
  {
    gmi_fail("EGADS model should only have one body");
  }

  *eg_body = eg_bodies[0];

  struct gmi_model *model;
  model = (struct gmi_model*)malloc(sizeof(*model));
  model->ops = &ops;

  EG_getBodyTopos(eg_body, NULL, NODE, &(model->n[0]), NULL);
  EG_getBodyTopos(eg_body, NULL, EDGE, &(model->n[1]), NULL);
  EG_getBodyTopos(eg_body, NULL, FACE, &(model->n[2]), NULL);
  EG_getBodyTopos(eg_body, NULL, SHELL, &(model->n[3]), NULL); // BODY?

  return model;
}
#else
static struct gmi_model* gmi_egads_load(const char* filename)
{
  (void)filename;
  /// TODO: chose a compile flag
  gmi_fail("recompile with -DUSE_EGADS=ON");
}
#endif

void gmi_egads_start(void)
{
  int status = EG_open(eg_context);
  if (status != EGADS_SUCCESS)
  {
    char str[50]; // big enough
    sprintf(str, "EGADS failed to open with error code: %d", status);
    gmi_fail(str);
  }
}

void gmi_egads_stop(void)
{
  EG_close(*eg_context);
}

void gmi_register_egads(void)
{
  ops.begin = begin;
  ops.next = next;
  ops.end = end;
  ops.dim = gmi_dim;
  ops.tag = gmi_tag;
  ops.find = find;
  ops.adjacent = adjacent;
  ops.eval = eval;
  ops.reparam = reparam;
  ops.periodic = periodic;
  ops.range = range;
  ops.closest_point = closest_point;
  ops.normal = normal;
  ops.first_derivative = first_derivative;
  ops.is_point_in_region = is_point_in_region;
  ops.is_in_closure_of = is_in_closure_of;
  ops.bbox = bbox;
  ops.is_discrete_ent = is_discrete_ent;
  ops.destroy = destroy;
  gmi_register(gmi_egads_load, "egads");
}
