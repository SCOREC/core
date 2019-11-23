/******************************************************************************

  Copyright 2014 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/
#include <gmi.h>
#include "egads.h"
#include "gmi_egads.h"
#include <math.h>

#include "gmi_egads_config.h"

// will be initialized by `gmi_egads_start`
ego eg_context;
// will be initialized by `gmi_egads_load`
ego eg_model;
// will be initialized by `gmi_egads_load`
ego eg_body;

struct egads_iter
{
   ego *ents;
   int nelem;
   int idx;
};

// struct gmi_iter* begin(struct gmi_model* m, int dim)
// {
//   printf("begin\n");
//   (void)m;
//   ego *eg_ents;
//   int nbodies = 0;
//   if (dim == 0)
//     EG_getBodyTopos(eg_body, NULL, NODE, &nbodies, &eg_ents);
//   else if (dim == 1)
//     EG_getBodyTopos(eg_body, NULL, EDGE, &nbodies, &eg_ents);
//   else if (dim == 2)
//     EG_getBodyTopos(eg_body, NULL, FACE, &nbodies, &eg_ents);
//   else if (dim == 3)
//     EG_getBodyTopos(eg_body, NULL, SHELL, &nbodies, &eg_ents); // BODY?
//   printf("got body topops of dim %d with %d bodies\n", dim, nbodies);
//   gmi_fail("oops\n");
//   if (dim >= 0 && dim <= 3)
//     return (struct gmi_iter*)eg_ents;
//   return (struct gmi_iter*)NULL;
// }

// struct gmi_ent* next(struct gmi_model* m, struct gmi_iter* i)
// {
//   printf("next\n");
//   (void)m;
//   ego *eg_ents = (ego*)i;
//   return (struct gmi_ent*)++eg_ents;
// }

// void end(struct gmi_model* m, struct gmi_iter* i)
// {
//   printf("end\n");
//   (void)m;
//   // I think this will create a memory leak as it won't free any of the
//   // values that came before
//   ego *eg_ents = (ego*)i;
//   EG_free(eg_ents);
// }

struct gmi_iter* begin(struct gmi_model* m, int dim)
{
  printf("begin\n");
  (void)m;
  ego *eg_ents;
  int nbodies = 0;
  if (dim == 0)
    EG_getBodyTopos(eg_body, NULL, NODE, &nbodies, &eg_ents);
  else if (dim == 1)
    EG_getBodyTopos(eg_body, NULL, EDGE, &nbodies, &eg_ents);
  else if (dim == 2)
    EG_getBodyTopos(eg_body, NULL, FACE, &nbodies, &eg_ents);
  else if (dim == 3)
    EG_getBodyTopos(eg_body, NULL, SHELL, &nbodies, &eg_ents); // BODY?
  printf("got body topops of dim %d with %d bodies\n", dim, nbodies);
  struct egads_iter *eg_iter;
  if (dim >= 0 && dim <= 3)
  {
    eg_iter = EG_alloc(sizeof(struct egads_iter));
    if (eg_iter == NULL)
    {
      printf("failed to make memory\n");
      gmi_fail("EG_alloc failed to allocate memory for iter");
      return NULL;
    }
    printf("malloc success\n");
    eg_iter->ents = eg_ents;
    eg_iter->nelem = nbodies;
    eg_iter->idx = 0;
    return (struct gmi_iter*)eg_iter;
  }
  return (struct gmi_iter*)NULL;
}

struct gmi_ent* next(struct gmi_model* m, struct gmi_iter* i)
{
  printf("next\n");
  (void)m;
  struct egads_iter *eg_iter = (struct egads_iter*)i;
  ego *eg_ent;
  if (eg_iter->idx < eg_iter->nelem)
  {
    eg_ent = &(eg_iter->ents[eg_iter->idx]);
    eg_iter->idx++;
    return (struct gmi_ent*)eg_ent;
  }
  else
    return NULL;
  // return (struct gmi_ent*)eg_ent;
}

void end(struct gmi_model* m, struct gmi_iter* i)
{
  printf("end\n");
  (void)m;
  // I think this will create a memory leak as it won't free any of the
  // values that came before
  // ego *eg_ents = (ego*)i;
  struct egads_iter *eg_iter = (struct egads_iter*)i;

  if (eg_iter != NULL)
  {
    EG_free(eg_iter->ents);
    EG_free(eg_iter);
  }
}

int get_dim(struct gmi_model* m, struct gmi_ent* e)
{
  printf("get dim\n");
  (void)m;
  ego *eg_ent = (ego*)e;
  if (eg_ent == NULL)
  {
    printf("null ptr\n");
  }
  int ent_type;
  int mtype;
  ego topref, prev, next;
  printf("above get info\n");
  int status = EG_getInfo(*eg_ent, &ent_type, &mtype, &topref, &prev, &next);
  printf("get info status: %d\n", status);
  printf("above get info\n");
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

int get_tag(struct gmi_model* m, struct gmi_ent* e)
{
  printf("tag\n");
  (void)m;
  ego *eg_ent = (ego*)e;
  return EG_indexBodyTopo(eg_body, *eg_ent);
}

struct gmi_ent* find(struct gmi_model* m, int dim, int tag)
{
  printf("find\n");
  (void)m;
  // ego *eg_ent = (ego*)e;
  ego eg_ent = NULL;
  if (dim == 0)
    EG_objectBodyTopo(eg_body, NODE, tag, &eg_ent);
  else if (dim == 1)
    EG_objectBodyTopo(eg_body, EDGE, tag, &eg_ent);
  else if (dim == 2)
    EG_objectBodyTopo(eg_body, FACE, tag, &eg_ent);
  else if (dim == 3)
    EG_objectBodyTopo(eg_body, SHELL, tag, &eg_ent);
  else
    gmi_fail("gmi_ent not found!");
  return (struct gmi_ent*)eg_ent;
}

struct gmi_set* adjacent(struct gmi_model* m, 
                                struct gmi_ent* e, 
                                int dim)
{
  printf("adjacent\n");
  (void)m;
  ego *eg_ent = (ego*)e;
  int num_adjacent = 0;
  ego *adjacent_ents = NULL;
  if (dim == 0)
    EG_getBodyTopos(eg_body, *eg_ent, NODE, &num_adjacent, &adjacent_ents);
  else if (dim == 1)
    EG_getBodyTopos(eg_body, *eg_ent, EDGE, &num_adjacent, &adjacent_ents);
  else if (dim == 2)
    EG_getBodyTopos(eg_body, *eg_ent, FACE, &num_adjacent, &adjacent_ents);
  else if (dim == 3)
    EG_getBodyTopos(eg_body, *eg_ent, SHELL, &num_adjacent, &adjacent_ents);
  
  struct gmi_set *gmi_adj_ent = gmi_make_set(num_adjacent);
  for (int i = 0; i < num_adjacent; ++i)
  {
    gmi_adj_ent->e[i] = (struct gmi_ent*)adjacent_ents[i];
  }
  EG_free(adjacent_ents);
  return gmi_adj_ent;
}

void eval(struct gmi_model* m, 
                 struct gmi_ent* e,
                 double const p[2],
                 double x[3])
{
  printf("eval\n");
  (void)m;
  double results[18];
  ego *eg_ent = (ego*)e;
  EG_evaluate(*eg_ent, p, results);
  x[0] = results[0];
  x[1] = results[1];
  x[2] = results[2];
}

void reparam(struct gmi_model* m,
                    struct gmi_ent* from,
                    double const from_p[2],
                    struct gmi_ent* to,
                    double to_p[2])
{
  printf("reparam\n");
  int from_dim, to_dim;
  from_dim = get_dim(m, from);
  to_dim = get_dim(m, to);
  ego *eg_from = (ego*)from;
  ego *eg_to = (ego*)to;
  if ((from_dim == 1) && (to_dim == 2))
  {
    EG_getEdgeUV(*eg_to, *eg_from, from_p[0], 1, to_p);
    return;
  }
  if ((from_dim == 0) && (to_dim == 2))
  {
    // Doesn't yet exist
    // EG_getVertexUV(*eg_to, *eg_from, to_p);
    gmi_fail("From node to surface reparam not implemented");
    return;
  }
  if ((from_dim == 0) && (to_dim == 1))
  {
    // Doesn't yet exist
    // EG_getVertexT(*eg_to, *eg_from, &to_p[0]);
    gmi_fail("From node to edge reparam not implemented");
    return;
  }
  gmi_fail("bad dimensions in gmi_egads reparam");
}

int periodic(struct gmi_model* m,
                    struct gmi_ent* e,
                    int dir)
{
  printf("periodic\n");
  int ent_dim = get_dim(m, e);
  int periodic;
  ego *eg_ent = (ego*)e;
  EG_getRange(*eg_ent, NULL, &periodic);

  if (dir == 1) // v direction
  {
    if (ent_dim == 2) // FACE
    {
      if (periodic == 0)
        return 0;
      if (periodic == 2)
        return 1;
    }
    else
      gmi_fail("v direction only exists for faces");
  }
  if (ent_dim == 1 || ent_dim == 2)
    return periodic;
  return 0;
}

void range(struct gmi_model* m,
                  struct gmi_ent* e,
                  int dir,
                  double r[2])
{
  printf("range\n");
  int ent_dim = get_dim(m, e);
  double range[4];
  ego *eg_ent = (ego*)e;
  EG_getRange(*eg_ent, range, NULL);
  if (dir == 1)
  {
    if (ent_dim == 2)
    {
      r[0] = range[2];
      r[1] = range[3];
    }
    else 
      gmi_fail("v direction only exists for faces");
  }
  else if (dir == 0)
  {
    r[0] = range[0];
    r[1] = range[1];
  }
}

void closest_point(struct gmi_model* m,
                          struct gmi_ent* e, 
                          double const from[3],
                          double to[3],
                          double to_p[2])
{
  printf("closest point\n");
  (void)m;
  ego *eg_ent = (ego*)e;
  double xyz[] = {from[0], from[1], from[2]};
  EG_invEvaluate(*eg_ent, &xyz[0], &to_p[0], &to[0]);
}

void normal(struct gmi_model* m,
                   struct gmi_ent* e,
                   double const p[2],
                   double n[3])
{
  printf("normal\n");
  double du[3], dv[3];
  m->ops->first_derivative(m, e, p, du, dv);
  // cross du and dv to get n
  n[0] = du[1]*dv[2] - du[2]*dv[1];
  n[1] = du[2]*dv[0] - du[0]*dv[2];
  n[2] = du[0]*dv[1] - du[1]*dv[0];

  int mtype = 0;
  ego *eg_ent = (ego*)e;
  // EG_getInfo(*eg_ent, NULL, &mtype, NULL, NULL, NULL);
  EG_getTopology(*eg_ent, NULL, NULL, &mtype, NULL, NULL, NULL, NULL);

  double n_mag = sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
  n[0] *= mtype / n_mag;
  n[1] *= mtype / n_mag;
  n[2] *= mtype / n_mag;
}

void first_derivative(struct gmi_model* m,
                             struct gmi_ent* e,
                             double const p[2],
                             double t0[3],
                             double t1[3])
{
  printf("first derivative\n");
  int ent_dim = get_dim(m, e);
  double results[18];
  ego *eg_ent = (ego*)e;
  EG_evaluate(*eg_ent, p, results);
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

int is_point_in_region(struct gmi_model* m,
                              struct gmi_ent* e,
                              double p[3])
{
  printf("is in region\n");
  (void)m;
  ego *eg_ent = (ego*)e;
  int status = EG_inTopology(*eg_ent, p);
  if (status == EGADS_SUCCESS)
    return 1;
  else
    return 0;
}

void bbox(struct gmi_model* m,
                 struct gmi_ent* e,
                 double bmin[3],
                 double bmax[3])
{
  printf("bbox\n");
  (void)m;
  double box[6];
  ego *eg_ent = (ego*)e;
  EG_getBoundingBox(*eg_ent, box);
  bmin[0] = box[0];
  bmin[1] = box[1];
  bmin[2] = box[2];
  bmax[0] = box[3];
  bmax[1] = box[4];
  bmax[2] = box[5];
}

/// For any given vertex, edge, or face, this function can be used
/// to see if the vertex/edge/face is adjacent to region.
int is_in_closure_of(struct gmi_model* m,
                            struct gmi_ent* e,
                            struct gmi_ent* et)
{
  printf("in closure of\n");
  ego *eg_ent = (ego*)e;
  ego *eg_region = (ego*)et;
  int ent_dim = get_dim(m, e);
  int num_adjacent = 0;
  ego *adjacent_ents = NULL;
  if (ent_dim == 0)
    EG_getBodyTopos(eg_body, *eg_region, NODE, &num_adjacent, &adjacent_ents);
  else if (ent_dim == 1)
    EG_getBodyTopos(eg_body, *eg_region, EDGE, &num_adjacent, &adjacent_ents);
  else if (ent_dim == 2)
    EG_getBodyTopos(eg_body, *eg_region, FACE, &num_adjacent, &adjacent_ents);
  else if (ent_dim == 3)
    EG_getBodyTopos(eg_body, *eg_region, SHELL, &num_adjacent, &adjacent_ents);
  for (int i = 0; i < num_adjacent; ++i)
  {
    if (EG_isEquivalent(*eg_ent, adjacent_ents[i]))
      return 1;
  }
  EG_free(adjacent_ents);
  return 0;
}

/// what does this function do?
int is_discrete_ent(struct gmi_model* m, struct gmi_ent* e)
{
  (void)m;
  (void)e;
  gmi_fail("is_discrete_ent not implemented");
}

void destroy(struct gmi_model* m)
{
  printf("destroy!\n");
  free(m);
}

struct gmi_model_ops ops;

/// TODO: Come up with a better flag? 
/// TODO: re-write for EGADSlite - model loading is different
// #ifdef HAVE_EGADS
#if 1
struct gmi_model* gmi_egads_load(const char* filename)
{
  printf("in gmi_egads_load\n");
  int load_status = EG_loadModel(eg_context, 0, filename, &eg_model);
  if (load_status != EGADS_SUCCESS)
  {
    char str[50]; // big enough
    sprintf(str, "EGADS failed to load model with error code: %d", load_status);
    gmi_fail(str);
  }
  printf("after EG_loadModel\n");

  /// TODO: only store the outputs I need, replace the rest with NULL
  int oclass, mtype, nbody, *senses;
  ego geom, *eg_bodies;
  int status = EG_getTopology(eg_model, &geom, &oclass, &mtype, NULL, &nbody,
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

  eg_body = eg_bodies[0];

  struct gmi_model *model;
  model = (struct gmi_model*)malloc(sizeof(*model));
  model->ops = &ops;

  EG_getBodyTopos(eg_body, NULL, NODE, &(model->n[0]), NULL);
  EG_getBodyTopos(eg_body, NULL, EDGE, &(model->n[1]), NULL);
  EG_getBodyTopos(eg_body, NULL, FACE, &(model->n[2]), NULL);
  // I believe this should be shell, but always seems to result in 1 shell
  EG_getBodyTopos(eg_body, NULL, SHELL, &(model->n[3]), NULL); // BODY?

  return model;
}
#else
struct gmi_model* gmi_egads_load(const char* filename)
{
  (void)filename;
  /// TODO: chose a compile flag
  gmi_fail("recompile with -DUSE_EGADS=ON");
}
#endif

void gmi_egads_start(void)
{
  printf("egads start\n");
  int status = EG_open(&eg_context);
  printf("after egads open\n");
  if (status != EGADS_SUCCESS)
  {
    char str[50]; // big enough
    sprintf(str, "EGADS failed to open with error code: %d", status);
    gmi_fail(str);
  }
}

void gmi_egads_stop(void)
{
  EG_close(eg_context);
}

void gmi_register_egads(void)
{
  ops.begin = begin;
  ops.next = next;
  ops.end = end;
  ops.dim = get_dim;
  ops.tag = get_tag;
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
