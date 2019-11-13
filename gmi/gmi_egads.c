/******************************************************************************

  Copyright 2014 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/
#include "gmi_egads.h"

// #include "egads.h"

// initialize to NULL, will be properly set by `gmi_egads_start`
ego *eg_context = NULL;
// initialize to NULL, will be properly set by `gmi_egads_load`
ego *eg_model = NULL;

static struct gmi_iter* begin(struct gmi_model* m, int dim)
{
  /// implement
}

static struct gmi_ent* next(struct gmi_model* m, struct gmi_iter* i)
{
  /// implement
}

static void end(struct gmi_model* m, struct gmi_iter* i)
{
  /// implement
}

static int dim(struct gmi_model* m, struct gmi_ent* e)
{
  /// implement
}

static int tag(struct gmi_model* m, struct gmi_ent* e)
{
  /// implement
}

static struct gmi_ent* find(struct gmi_model* m, int dim, int tag)
{
  /// implement
}

static struct gmi_set* adjacent(struct gmi_model* m, 
                                struct gmi_ent* e, 
                                int dim)
{
  /// implement
}

static void eval(struct gmi_model* m, 
                 struct gmi_ent* e,
                 double const p[2],
                 double x[3])
{
  /// implement
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
  /// implement
}

static void range(struct gmi_model* m,
                  struct gmi_ent* e,
                  int dim,
                  double r[2])
{
  /// implement
}

static void closest_point(struct gmi_model* m,
                          struct gmi_ent* e, 
                          double const from[3],
                          double to[3],
                          double to_p[2])
{
  /// implement
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
  /// implement
}

static int is_point_in_region(struct gmi_model* m,
                              struct gmi_ent* e,
                              double p[3])
{
  /// implement
}

static void bbox(struct gmi_model* m,
                 struct gmi_ent* e,
                 double bmin[3],
                 double bmax[3])
{
  /// implement
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

  struct gmi_model *model;
  model = (struct gmi_model*)malloc(sizeof(*model));
  model->ops = &ops;

  model->n[0] = EG_getBodyNumNodes(eg_model); // function call for num vertices on geo model
  model->n[1] = EG_getBodyNumEdges(eg_model); // function call for num edges on geo model
  model->n[2] = EG_getBodyNumFaces(eg_model); // function call for num faces on geo model
  model->n[3] = EG_getBodyNumShells(eg_model); // function call for num regions on geo model

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
  ops.dim = dim;
  ops.tag = tag;
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
