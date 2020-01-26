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

struct egads_ent
{
  ego *ego_ent;
  int dim;
  int tag;
} typedef egads_ent;

struct egads_iter
{
   egads_ent *ents;
   int nelem;
   int idx;
};

int ****adjacency_table;


/// TODO: consider optimizing adjacency tables and access
void get_3D_adjacency(struct gmi_model* m,
                      egads_ent* ent, 
                      int adj_dim, 
                      int *num_adjacent, 
                      egads_ent** adjacent_ents)
{
  int ent_dim = m->ops->dim(m, (struct gmi_ent*)ent);
  int ent_tag = m->ops->tag(m, (struct gmi_ent*)ent);
  
  int *adj_tags = adjacency_table[adj_dim][ent_dim][ent_tag];
  *num_adjacent = adj_tags[0];
  *adjacent_ents = (egads_ent*)EG_calloc(*num_adjacent, sizeof(egads_ent*));

  if (adj_dim == 3)
  {
    for (int i = 0; i < *num_adjacent; i++)
    {
      adjacent_ents[i]->ego_ent = NULL;
      adjacent_ents[i]->dim = 3;
      adjacent_ents[i]->tag = adj_tags[i+1]; // first entry is the number of adjacent
    }
  }
  else
  {
    for (int i = 0; i < *num_adjacent; i++)
    {
      egads_ent *eg_ent = (egads_ent*)m->ops->find(m, adj_dim, adj_tags[i]);
      adjacent_ents[i]->ego_ent = eg_ent->ego_ent;
      adjacent_ents[i]->dim = -1;
      adjacent_ents[i]->tag = -1;
    }
  }
}

/// TODO: implement based on adjacent face's bounding boxes
void get_3D_bounding_box(egads_ent *ent, double *box)
{
  (void)ent;
  (void)box;
}

/// reparameterize a vertex onto an edge
void getVertexT(struct gmi_model* m, struct gmi_ent* to, struct gmi_ent* from, double* t)
{
  double diff;
  double t_range[2];
  m->ops->range(m, to, 0, &(t_range[0]));
  printf("got range\n");
  double vtx_pnt[3];
  double p[] = {0, 0};
  m->ops->eval(m, from, p, &(vtx_pnt[0]));
  printf("eval 1\n");
  double t_pnt[3];
  m->ops->eval(m, to, &(t_range[0]), &(t_pnt[0]));
  printf("eval 2\n");
  diff = sqrt(pow(vtx_pnt[0] - t_pnt[0], 2) + 
              pow(vtx_pnt[1] - t_pnt[1], 2) + 
              pow(vtx_pnt[2] - t_pnt[2], 2));
  printf("diff 1: %f\n", diff);
  if (diff < 0.001)
  {
    *t = t_range[0];
  }
  else
  {
    m->ops->eval(m, to, &(t_range[1]), &(t_pnt[0]));
    diff = sqrt(pow(vtx_pnt[0] - t_pnt[0], 2) + 
                pow(vtx_pnt[1] - t_pnt[1], 2) + 
                pow(vtx_pnt[2] - t_pnt[2], 2));
    printf("diff (if here should be small): %f", diff);
    *t = t_range[1];
  }
  return;
}

// void getVertexUV(const ego to, const ego from, double to_p[2])
// {

// }

struct gmi_iter* begin(struct gmi_model* m, int dim)
{
  printf("begin\n");
  (void)m;
  ego *ego_ents;
  int nbodies = m->n[3];
  if (dim == 0)
    EG_getBodyTopos(eg_body, NULL, NODE, &nbodies, &ego_ents);
  else if (dim == 1)
    EG_getBodyTopos(eg_body, NULL, EDGE, &nbodies, &ego_ents);
  else if (dim == 2)
    EG_getBodyTopos(eg_body, NULL, FACE, &nbodies, &ego_ents);

  egads_ent *eg_ents = (egads_ent*)EG_alloc(nbodies*sizeof(egads_ent*));
  for (int i = 0; i < nbodies; i++)
  {
    if (dim == 3)
    {
      eg_ents[i].ego_ent = NULL;
      eg_ents[i].dim = 3;
      eg_ents[i].tag = i;
    }
    else
    {
      eg_ents[i].ego_ent = &ego_ents[i];
      eg_ents[i].dim = -1;
      eg_ents[i].tag = -1;
    }
  }
  struct egads_iter *eg_iter;
  if (dim >= 0 && dim <= 3)
  {
    eg_iter = EG_alloc(sizeof(struct egads_iter));
    if (eg_iter == NULL)
    {
      gmi_fail("EG_alloc failed to allocate memory for iter");
      return NULL;
    }
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
  // ego *eg_ent;
  if (eg_iter->idx < eg_iter->nelem)
  {
    // eg_ent = &(eg_iter->ents[eg_iter->idx]);
    // eg_iter->idx++;
    // ego ent = *eg_ent;
    // printf("dereferneced!\n");
    // printf("ego oclass: %d\n", ent->oclass);
    // return (struct gmi_ent*)eg_ent;
    return (struct gmi_ent*)&(eg_iter->ents[eg_iter->idx++]);
  }
  else
    return NULL;
  // return (struct gmi_ent*)eg_ent;
}

void end(struct gmi_model* m, struct gmi_iter* i)
{
  printf("end\n");
  (void)m;
  (void)i;
  // I think this will create a memory leak as it won't free any of the
  // values that came before
  // ego *eg_ents = (ego*)i;
  // struct egads_iter *eg_iter = (struct egads_iter*)i;

  /// I think freeing the array here will free it too early, 
  // if (eg_iter != NULL)
  // {
  //   EG_free(eg_iter->ents);
  //   EG_free(eg_iter);
  // }
}

int get_dim(struct gmi_model* m, struct gmi_ent* e)
{
  printf("get dim\n");
  (void)m;
  egads_ent* eg_ent = (egads_ent*)e;
  if (eg_ent->dim != -1) // 3D entity
  {
    return eg_ent->dim;
  }
  ego *ego_ent = eg_ent->ego_ent;
  if (ego_ent == NULL)
  {
    printf("null ptr\n");
  }
  // printf("derefence....\n");
  // ego ent = *ego_ent;
  // printf("dereferneced!\n");
  // printf("ego oclass: %d\n", ent->oclass);
  int oclass;
  int mtype;
  ego topref, prev, next;
  // printf("above get info\n");
  int status = EG_getInfo(*ego_ent, &oclass, &mtype, &topref, &prev, &next);

  if (status != EGADS_SUCCESS)
  {
    printf("error!\n");
  }
  // printf("get info status: %d\n", status);

  if (oclass == NODE)
    return 0;
  else if (oclass == EDGE)
    return 1;
  else if (oclass == FACE)
    return 2;
  return -1;
}

int get_tag(struct gmi_model* m, struct gmi_ent* e)
{
  printf("tag\n");
  (void)m;
  egads_ent* eg_ent = (egads_ent*)e;
  if (eg_ent->tag != -1) // 3D entity
  {
    return eg_ent->tag;
  }
  ego *ego_ent = eg_ent->ego_ent;
  return EG_indexBodyTopo(eg_body, *ego_ent);
}

struct gmi_ent* find(struct gmi_model* m, int dim, int tag)
{
  printf("find\n");
  (void)m;

  /// Not sure if this is the best way to handle this, previously was returning
  /// address to stack memory, so when memory dereferenced was not an ego
  /// might need to think about when to free this memory.
  egads_ent *eg_ent = (egads_ent*)EG_alloc(sizeof(egads_ent*));
  eg_ent->dim = -1;
  eg_ent->tag = -1;

  if (dim == 0)
    EG_objectBodyTopo(eg_body, NODE, tag, eg_ent->ego_ent);
  else if (dim == 1)
    EG_objectBodyTopo(eg_body, EDGE, tag, eg_ent->ego_ent);
  else if (dim == 2)
    EG_objectBodyTopo(eg_body, FACE, tag, eg_ent->ego_ent);
  else if (dim == 3)
  {
    eg_ent->dim = dim;
    eg_ent->tag = tag;
  }
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
  egads_ent *eg_ent = (egads_ent*)e;
  ego *ego_ent = eg_ent->ego_ent;
  int num_adjacent = 0;
  egads_ent *adjacent_ents;
  ego *adjacent_egos;
  if (eg_ent->dim == 3 || dim == 3)
  {
    get_3D_adjacency(m, eg_ent, dim, &num_adjacent, &adjacent_ents);
  }
  else // only dealing with egos
  {
    if (dim == 0)
      EG_getBodyTopos(eg_body, *ego_ent, NODE, &num_adjacent, &adjacent_egos);
    else if (dim == 1)
      EG_getBodyTopos(eg_body, *ego_ent, EDGE, &num_adjacent, &adjacent_egos);
    else if (dim == 2)
      EG_getBodyTopos(eg_body, *ego_ent, FACE, &num_adjacent, &adjacent_egos);

    adjacent_ents = (egads_ent*)EG_calloc(num_adjacent, sizeof(egads_ent*));
    for (int i = 0; i < num_adjacent; i++)
    {
      adjacent_ents[i].ego_ent = &adjacent_egos[i];
      adjacent_ents[i].dim = -1;
      adjacent_ents[i].tag = -1;
    }
  }

  struct gmi_set *gmi_adj_ent = gmi_make_set(num_adjacent);
  for (int i = 0; i < num_adjacent; ++i)
  {
    egads_ent* adj_ent = &(adjacent_ents[i]);
    gmi_adj_ent->e[i] = (struct gmi_ent*)adj_ent;
    // *(gmi_adj_ent->e + i) = (struct gmi_ent*)adjacent_ents[i];
    // printf("adjacent ent: %d oclass: %d\n", i, adjacent_ents[i]->oclass);
  }
  // EG_free(adjacent_ents);

  // // maybe by doing [] I loose the pointer type, need to think carefully about it
  // // struct gmi_ent* gent = (struct gmi_ent*)gmi_adj_ent->e[0];
  // ego* adj_ent = &(adjacent_ents[0]);
  // struct gmi_ent* gent = (struct gmi_ent*)adj_ent;
  // printf("cast 1\n");
  // ego* eg_ent2 = (ego*)gent;
  // printf("cast 2\n");
  // ego eg_ent3 = *eg_ent2;
  // printf("cast 3\n");
  // printf("oclass round 2: %d\n", eg_ent3->oclass);
  return gmi_adj_ent;
}

/// TODO: error check to make sure ego_ent != NULL?
void eval(struct gmi_model* m, 
          struct gmi_ent* e,
          double const p[2],
          double x[3])
{
  printf("eval\n");
  // (void)m;
  double results[18];
  egads_ent *eg_ent = (egads_ent*)e;
  ego *ego_ent = eg_ent->ego_ent;
  int dim = m->ops->dim(m, e);
  // printf("dim: %d\n", dim);
  if (dim > 0)
  {
    EG_evaluate(*ego_ent, p, results);
    x[0] = results[0];
    x[1] = results[1];
    x[2] = results[2];
  }
  else if (dim == 0)
  {
    // int mtype;
    double data[4];
    int oclass, mtype, nbody, *senses;
    ego geom, *eg_bodies;
    // EG_getTopology(eg_model, &geom, &oclass, &mtype, NULL, &nbody,
    //                       &eg_bodies, &senses);
    EG_getTopology(*ego_ent, &geom, &oclass, &mtype, data, &nbody, &eg_bodies, &senses);
    // printf("after get topo\n");
    x[0] = data[0];
    x[1] = data[1];
    x[2] = data[2];
  }
}

/// TODO: error check to make sure ego_ent != NULL?
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
  egads_ent *eg_ent_from = (egads_ent*)from;
  ego *ego_from = eg_ent_from->ego_ent;

  egads_ent *eg_ent_to = (egads_ent*)to;
  ego *ego_to = eg_ent_to->ego_ent;

  if ((from_dim == 1) && (to_dim == 2))
  {
    EG_getEdgeUV(*ego_to, *ego_from, 1, from_p[0], to_p);
    return;
  }
  if ((from_dim == 0) && (to_dim == 2))
  {
    printf("reparam from %d to %d not implemented\n", from_dim, to_dim);
    // Doesn't yet exist
    // EG_getVertexUV(*ego_to, *ego_from, to_p);
    gmi_fail("From node to surface reparam not implemented");
    return;
  }
  if ((from_dim == 0) && (to_dim == 1))
  {
    printf("reparam from %d to %d not implemented\n", from_dim, to_dim);
    // Doesn't yet exist
    getVertexT(m, to, from, &to_p[0]);
    // gmi_fail("From node to edge reparam not implemented");
    return;
  }
  printf("attempted reparam from %d to %d\n", from_dim, to_dim);
  gmi_fail("bad dimensions in gmi_egads reparam");
}

/// TODO: error check to make sure ego_ent != NULL?
int periodic(struct gmi_model* m,
             struct gmi_ent* e,
             int dir)
{
  printf("periodic\n");
  int ent_dim = get_dim(m, e);
  int periodic;
  egads_ent *eg_ent = (egads_ent*)e;
  ego *ego_ent = eg_ent->ego_ent;

  double range[4];
  EG_getRange(*ego_ent, range, &periodic);

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

/// TODO: error check to make sure ego_ent != NULL?
void range(struct gmi_model* m,
           struct gmi_ent* e,
           int dir,
           double r[2])
{
  printf("range\n");
  int ent_dim = get_dim(m, e);
  double range[4];
  int periodic;
  egads_ent *eg_ent = (egads_ent*)e;
  ego *ego_ent = eg_ent->ego_ent;

  EG_getRange(*ego_ent, range, &periodic);
  printf("after EG_getRange\n");
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

/// TODO: error check to make sure ego_ent != NULL?
void closest_point(struct gmi_model* m,
                   struct gmi_ent* e, 
                   double const from[3],
                   double to[3],
                   double to_p[2])
{
  printf("closest point\n");
  (void)m;
  egads_ent *eg_ent = (egads_ent*)e;
  ego *ego_ent = eg_ent->ego_ent;
  double xyz[] = {from[0], from[1], from[2]};
  EG_invEvaluate(*ego_ent, &xyz[0], &to_p[0], &to[0]);
}

/// TODO: error check to make sure ego_ent != NULL?
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

  // int mtype = 0;
  egads_ent *eg_ent = (egads_ent*)e;
  ego *ego_ent = eg_ent->ego_ent;
  // EG_getInfo(*ego_ent, NULL, &mtype, NULL, NULL, NULL);
  // EG_getTopology(*ego_ent, NULL, NULL, &mtype, NULL, NULL, NULL, NULL);
  double data[4];
  int oclass, mtype, nbody, *senses;
  ego geom, *eg_bodies;
  // EG_getTopology(eg_model, &geom, &oclass, &mtype, NULL, &nbody,
  //                       &eg_bodies, &senses);
  EG_getTopology(*ego_ent, &geom, &oclass, &mtype, data, &nbody, &eg_bodies, &senses);

  double n_mag = sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
  n[0] *= mtype / n_mag;
  n[1] *= mtype / n_mag;
  n[2] *= mtype / n_mag;
}

/// TODO: error check to make sure ego_ent != NULL?
void first_derivative(struct gmi_model* m,
                      struct gmi_ent* e,
                      double const p[2],
                      double t0[3],
                      double t1[3])
{
  printf("first derivative\n");
  int ent_dim = get_dim(m, e);
  double results[18];
  egads_ent *eg_ent = (egads_ent*)e;
  ego *ego_ent = eg_ent->ego_ent;
  EG_evaluate(*ego_ent, p, results);
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

/// TODO: make this work for new 3D object
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

/// TODO: make this work for new 3D object
void bbox(struct gmi_model* m,
          struct gmi_ent* e,
          double bmin[3],
          double bmax[3])
{
  printf("bbox\n");
  (void)m;
  double box[6];
  egads_ent *eg_ent = (egads_ent*)e;
  ego *ego_ent = eg_ent->ego_ent;

  if (eg_ent->dim != -1)
  {
    get_3D_bounding_box(eg_ent, box);
  }
  else
  {
    EG_getBoundingBox(*ego_ent, box);
  }
  bmin[0] = box[0];
  bmin[1] = box[1];
  bmin[2] = box[2];
  bmax[0] = box[3];
  bmax[1] = box[4];
  bmax[2] = box[5];
}

/// For any given vertex, edge, or face e, this function can be used
/// to see if the e is in the closure of entity et.
int is_in_closure_of(struct gmi_model* m,
                     struct gmi_ent* e,
                     struct gmi_ent* et)
{
  printf("in closure of\n");
  egads_ent *eg_ent = (egads_ent*)e;
  ego *ego_ent = eg_ent->ego_ent;

  egads_ent *eg_region = (egads_ent*)et;
  ego *ego_region = eg_ent->ego_ent;

  int num_adjacent = 0;

  int ent_dim = get_dim(m, e);
  int region_dim = get_dim(m, et);
  if (ent_dim == 3 || region_dim == 3)
  {
    egads_ent *adjacent_ents;
    // get entities of dim ent_dim adjacent to eg_region (will be downward adjacent)
    get_3D_adjacency(m, eg_region, ent_dim, &num_adjacent, &adjacent_ents);

    for (int i = 0; i < num_adjacent; i++)
    {
      if (adjacent_ents[i].tag == eg_ent->tag)
      {
        EG_free(adjacent_ents);
        return 1;
      }
    }
    EG_free(adjacent_ents);
    return 0;
  }
  else
  {
    ego *adjacent_egos = NULL;
    if (ent_dim == 0)
      EG_getBodyTopos(eg_body, *ego_region, NODE, &num_adjacent, &adjacent_egos);
    else if (ent_dim == 1)
      EG_getBodyTopos(eg_body, *ego_region, EDGE, &num_adjacent, &adjacent_egos);
    else if (ent_dim == 2)
      EG_getBodyTopos(eg_body, *ego_region, FACE, &num_adjacent, &adjacent_egos);
    for (int i = 0; i < num_adjacent; ++i)
    {
      if (EG_isEquivalent(*ego_ent, adjacent_egos[i]))
      {
        EG_free(adjacent_egos);
        return 1;
      }
    }
    EG_free(adjacent_egos);
    return 0;
  }
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
