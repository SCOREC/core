/******************************************************************************

  Copyright 2014 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/
#include <math.h>
#include <string.h>
#include <stdio.h>

#include "egads.h"

#include <gmi.h>
#include "gmi_egads.h"
#include "gmi_egads_config.h"

// will be initialized by `gmi_egads_start`
ego eg_context;
// will be initialized by `gmi_egads_load`
ego eg_model;
// will be initialized by `gmi_egads_load`
ego eg_body;

typedef struct egads_ent
{
  ego ego_ent;
  int dim;
  int tag;
} egads_ent;

/// global array of model entities
/// egads_global_ents[0] - array of model verts
/// egads_global_ents[1] - array of model edges
/// egads_global_ents[2] - array of model faces
/// egads_global_ents[3] - array of model regions
egads_ent *egads_global_ents[4];

struct egads_iter
{
   egads_ent *ents;
   int nelem;
   int idx;
};

/// adjacency table to be used for the "dummy" 3D elements that EGADS doesn't
/// natively track
int **adjacency_table[6];

/// read adjacency table from binary file
void read_adj_table(const char* filename,
                    int *nregions)
{
  char *sup_filename;
  sup_filename = EG_alloc(strlen(filename)+4+1);
  if (sup_filename == NULL)
  {
    gmi_fail("failed to allocate memory for new string");
  }
  sup_filename[0] = '\0';
  strcat(sup_filename,filename);
  strcat(sup_filename, ".sup");

  FILE *adj_table_file = fopen(sup_filename, "rb");
  if (adj_table_file == NULL)
  {
    gmi_fail("failed to open supplementary EGADS model file!");
  }

  int header[6];
  fread(header, sizeof(int), 6, adj_table_file);

  *nregions = header[0];

  for (int i = 0; i < 6; ++i)
  {
    adjacency_table[i] = (int**)EG_alloc(sizeof(*(adjacency_table[i]))
                                         * header[i]);
    if (adjacency_table[i] == NULL) {
      char fail[50];
      sprintf(fail, "failed to alloc memory for adjacency_table[%d]", i);
      /// TODO: this could cause memory leak
      gmi_fail(fail);
    }
    for (int j = 0; j < header[i]; ++j) {
      int nadjacent = -1;
      fread(&nadjacent, sizeof(int), 1, adj_table_file);
      adjacency_table[i][j] = (int*)EG_alloc(sizeof(*(adjacency_table[i][j]))
                                             * (nadjacent+1));
      if (adjacency_table[i][j] == NULL) {
        char fail[50];
        sprintf(fail, "failed to alloc memory for "
                "adjacency_table[%d][%d]", i,j);
        /// TODO: this could cause memory leak
        gmi_fail(fail);
      }
      adjacency_table[i][j][0] = nadjacent;
      fread(&(adjacency_table[i][j][1]), sizeof(int), nadjacent,
              adj_table_file);
    }
  }
  fclose(adj_table_file);
}


/// TODO: consider optimizing adjacency tables and access
void get_3D_adjacency(struct gmi_model* m,
                      egads_ent* ent, 
                      int adj_dim, 
                      int *num_adjacent, 
                      egads_ent*** adjacent_ents)
{
  int ent_dim = m->ops->dim(m, (struct gmi_ent*)ent);
  int ent_tag = m->ops->tag(m, (struct gmi_ent*)ent);
  
  int pairing = -1;
  if (adj_dim == 0 && ent_dim == 3)
    pairing = 0;
  else if (adj_dim == 1 && ent_dim == 3)
    pairing = 1;
  else if (adj_dim == 2 && ent_dim == 3)
    pairing = 2;
  else if (adj_dim == 3 && ent_dim == 0)
    pairing = 3;
  else if (adj_dim == 3 && ent_dim == 1)
    pairing = 4;
  else if (adj_dim == 3 && ent_dim == 2)
    pairing = 5;
  else
    gmi_fail("bad dims in get_3D_adjacency!");


  int *adj_tags = adjacency_table[pairing][ent_tag-1];
  *num_adjacent = adj_tags[0];
  *adjacent_ents = (egads_ent**)EG_alloc(sizeof(**adjacent_ents) * (*num_adjacent));

  for (int i = 0; i < *num_adjacent; ++i)
  {
    (*adjacent_ents)[i] = (egads_ent*)m->ops->find(m, adj_dim, adj_tags[i+1]);
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
    printf("diff (if here should be small): %f\n", diff);
    *t = t_range[1];
  }
  return;
}

/// TODO: handle degenerate edges (cant eval on them)
/// function to reparameterize a vertex onto a face
/// (get the u,v parametric coodinates associated with the vertex)
void getVertexUV(struct gmi_model* m, struct gmi_ent* to,
                 struct gmi_ent* from, double to_p[2])
{
  struct gmi_set* adj_faces;
  struct gmi_set* adj_edges = gmi_adjacent(m, from, 1);

  for (int i = 0; i < adj_edges->n; ++i)
  {
    adj_faces = gmi_adjacent(m, adj_edges->e[i], 2);
    for (int j = 0; j < adj_faces->n; ++j)
    {
      if (adj_faces->e[j] == to)
      {
        double t;
        getVertexT(m, adj_edges->e[i], from, &t);
        m->ops->reparam(m, adj_edges->e[i], &t, to, to_p);
        goto cleanup;
      }
    }
    gmi_free_set(adj_faces);
  }

  cleanup:
    gmi_free_set(adj_faces);
    gmi_free_set(adj_edges);
    return;
}

struct gmi_iter* begin(struct gmi_model* m, int dim)
{
  printf("begin\n");
  
  struct egads_iter *eg_iter;
  if (dim >= 0 && dim <= 3)
  {
    eg_iter = EG_alloc(sizeof(*eg_iter));
    if (eg_iter == NULL)
    {
      gmi_fail("EG_alloc failed to allocate memory for iter");
      return (struct gmi_iter*)NULL;
    }
    int nents = m->n[dim];
    eg_iter->ents = &(egads_global_ents[dim][0]);
    eg_iter->nelem = nents;
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
  if (eg_iter->idx < eg_iter->nelem)
  {
    return (struct gmi_ent*)(eg_iter->ents+eg_iter->idx++);
  }
  return (struct gmi_ent*)NULL;
}

void end(struct gmi_model* m, struct gmi_iter* i)
{
  printf("end\n");
  (void)m;

  struct egads_iter *eg_iter = (struct egads_iter*)i;
  /// I think freeing the array here will free it too early, 
  if (eg_iter != NULL)
  {
    EG_free(eg_iter);
  }
}

int get_dim(struct gmi_model* m, struct gmi_ent* e)
{
  printf("get dim\n");
  (void)m;
  egads_ent* eg_ent = (egads_ent*)e;
  return eg_ent->dim;
}

int get_tag(struct gmi_model* m, struct gmi_ent* e)
{
  printf("get tag\n");
  (void)m;
  egads_ent* eg_ent = (egads_ent*)e;
  return eg_ent->tag;
}

struct gmi_ent* find(struct gmi_model* m, int dim, int tag)
{
  printf("find\n");
  (void)m;
  return (struct gmi_ent*)&(egads_global_ents[dim][tag-1]);
}

// find (dim)-dimensional entities adjacent to e
struct gmi_set* adjacent(struct gmi_model* m, 
                         struct gmi_ent* e, 
                         int dim)
{
  egads_ent *eg_ent = (egads_ent*)e;
  printf("adjacent! dim: %d, ent dim: %d, ent tag: %d\n", dim, eg_ent->dim, eg_ent->tag);

  int num_adjacent = 0;

  egads_ent **adjacent_ents = NULL;
  ego *adjacent_egos;

  if (eg_ent->dim == 3 || dim == 3)
  {
    get_3D_adjacency(m, eg_ent, dim, &num_adjacent, &adjacent_ents);
  }
  else // only dealing with egos
  {
    if (dim == 0)
      EG_getBodyTopos(eg_body, (eg_ent->ego_ent), 20, &num_adjacent, &adjacent_egos);
    else if (dim == 1)
      EG_getBodyTopos(eg_body, (eg_ent->ego_ent), 21, &num_adjacent, &adjacent_egos);
    else if (dim == 2)
      EG_getBodyTopos(eg_body, (eg_ent->ego_ent), 23, &num_adjacent, &adjacent_egos);

    adjacent_ents = (egads_ent**)EG_alloc(sizeof(*adjacent_ents) * num_adjacent);
    for (int i = 0; i < num_adjacent; i++)
    {
      int adj_ent_tag = EG_indexBodyTopo(eg_body, adjacent_egos[i]);
      // adjacent_ents[i].ego_ent = adjacent_egos[i];
      // adjacent_ents[i].dim = dim;
      // adjacent_ents[i].tag = adj_ent_tag;
      adjacent_ents[i] = (egads_ent*)m->ops->find(m, dim, adj_ent_tag);
    }
  }

  struct gmi_set *gmi_adj_ent = gmi_make_set(num_adjacent);
  for (int i = 0; i < num_adjacent; ++i)
  {
    // int tag = m->ops->tag(m, (struct gmi_ent*)&(adjacent_ents[i]));
    gmi_adj_ent->e[i] = (struct gmi_ent*)adjacent_ents[i];
  }
  EG_free(adjacent_ents);
  adjacent_ents = NULL;

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
  ego ego_ent = eg_ent->ego_ent;
  int dim = m->ops->dim(m, e);
  // printf("dim: %d\n", dim);
  if (dim > 0 && dim < 3)
  {
    EG_evaluate(ego_ent, p, results);
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
    EG_getTopology(ego_ent, &geom, &oclass, &mtype, data, &nbody,
                   &eg_bodies, &senses);
    // printf("after get topo\n");
    x[0] = data[0];
    x[1] = data[1];
    x[2] = data[2];
  }
  else if (dim == 3)
  {
    gmi_fail("cannot eval 3D entity!");
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
  ego ego_from = eg_ent_from->ego_ent;

  egads_ent *eg_ent_to = (egads_ent*)to;
  ego ego_to = eg_ent_to->ego_ent;

  if ((from_dim == 1) && (to_dim == 2))
  {
    printf("reparam from %d to %d\n", from_dim, to_dim);
    EG_getEdgeUV(ego_to, ego_from, 1, from_p[0], to_p);
    return;
  }
  if ((from_dim == 0) && (to_dim == 2))
  {
    printf("reparam from %d to %d\n", from_dim, to_dim);
    // getVertexUV(*ego_to, *ego_from, to_p);
    getVertexUV(m, to, from, to_p);
    // gmi_fail("From node to surface reparam not implemented");
    return;
  }
  if ((from_dim == 0) && (to_dim == 1))
  {
    printf("reparam from %d to %d\n", from_dim, to_dim);
    getVertexT(m, to, from, &to_p[0]);
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
  ego ego_ent = eg_ent->ego_ent;

  double range[4];
  EG_getRange(ego_ent, range, &periodic);

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
  ego ego_ent = eg_ent->ego_ent;

  EG_getRange(ego_ent, range, &periodic);
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
  ego ego_ent = eg_ent->ego_ent;
  double xyz[] = {from[0], from[1], from[2]};
  EG_invEvaluate(ego_ent, &xyz[0], &to_p[0], &to[0]);
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
  ego ego_ent = eg_ent->ego_ent;
  // EG_getInfo(*ego_ent, NULL, &mtype, NULL, NULL, NULL);
  // EG_getTopology(*ego_ent, NULL, NULL, &mtype, NULL, NULL, NULL, NULL);
  double data[4];
  int oclass, mtype, nbody, *senses;
  ego geom, *eg_bodies;
  // EG_getTopology(eg_model, &geom, &oclass, &mtype, NULL, &nbody,
  //                       &eg_bodies, &senses);
  EG_getTopology(ego_ent, &geom, &oclass, &mtype, data, &nbody, &eg_bodies, &senses);

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
  ego ego_ent = eg_ent->ego_ent;
  EG_evaluate(ego_ent, p, results);
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
  egads_ent *eg_ent = (egads_ent*)e;
  if (eg_ent->ego_ent == NULL)
  {
    /// TODO: implement
    return -1;
  }
  else
  {
    ego ego_ent = eg_ent->ego_ent;
    int status = EG_inTopology(ego_ent, p);
    if (status == EGADS_SUCCESS)
      return 1;
    else
      return 0;
  }
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
  ego ego_ent = eg_ent->ego_ent;

  if (eg_ent->ego_ent == NULL)
  {
    get_3D_bounding_box(eg_ent, box);
  }
  else
  {
    EG_getBoundingBox(ego_ent, box);
  }
  bmin[0] = box[0];
  bmin[1] = box[1];
  bmin[2] = box[2];
  bmax[0] = box[3];
  bmax[1] = box[4];
  bmax[2] = box[5];
}


/// TODO: seems like this should call adjacent?
/// For any given vertex, edge, or face e, this function can be used
/// to see if the e is in the closure of entity et.
int is_in_closure_of(struct gmi_model* m,
                     struct gmi_ent* e,
                     struct gmi_ent* et)
{
  printf("in closure of\n");
  egads_ent *eg_ent = (egads_ent*)e;
  ego ego_ent = eg_ent->ego_ent;

  egads_ent *eg_region = (egads_ent*)et;
  ego ego_region = eg_ent->ego_ent;

  int num_adjacent = 0;

  int ent_dim = get_dim(m, e);
  int region_dim = get_dim(m, et);
  if (ent_dim == 3 || region_dim == 3)
  {
    egads_ent **adjacent_ents;
    // get entities of dim ent_dim adjacent to eg_region (will be downward adjacent)
    get_3D_adjacency(m, eg_region, ent_dim, &num_adjacent, &adjacent_ents);

    for (int i = 0; i < num_adjacent; i++)
    {
      if (adjacent_ents[i]->tag == eg_ent->tag)
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
      EG_getBodyTopos(eg_body, ego_region, NODE, &num_adjacent, &adjacent_egos);
    else if (ent_dim == 1)
      EG_getBodyTopos(eg_body, ego_region, EDGE, &num_adjacent, &adjacent_egos);
    else if (ent_dim == 2)
      EG_getBodyTopos(eg_body, ego_region, FACE, &num_adjacent, &adjacent_egos);
    for (int i = 0; i < num_adjacent; ++i)
    {
      if (EG_isEquivalent(ego_ent, adjacent_egos[i]))
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

/// TODO: free adjacency table too
void destroy(struct gmi_model* m)
{
  printf("destroy!\n");
  for (int i = 0; i < 4; ++i)
  {
    EG_free(egads_global_ents[i]);
  }

  int sizes[] = {m->n[3], m->n[3], m->n[3],
                 m->n[0], m->n[1], m->n[2]};
  for (int i = 0; i < 6; ++i)
  {
    for (int j = 0; j < sizes[i]; ++j)
    {
      EG_free(adjacency_table[i][j]);
    }
    EG_free(adjacency_table[i]);
  }

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

  int nregions = 1; // read adjacency file to find this number
  read_adj_table(filename, &nregions);
  // eg_body = eg_bodies[0];

  // struct gmi_model *model;
  // model = (struct gmi_model*)malloc(sizeof(*model));
  // model->ops = &ops;

  // EG_getBodyTopos(eg_body, NULL, NODE, &(model->n[0]), NULL);
  // EG_getBodyTopos(eg_body, NULL, EDGE, &(model->n[1]), NULL);
  // EG_getBodyTopos(eg_body, NULL, FACE, &(model->n[2]), NULL);
  // // I believe this should be shell, but always seems to result in 1 shell
  // EG_getBodyTopos(eg_body, NULL, SHELL, &(model->n[3]), NULL); // BODY?

  return gmi_egads_init(eg_bodies[0], nregions);
}
#else
struct gmi_model* gmi_egads_load(const char* filename)
{
  (void)filename;
  /// TODO: chose a compile flag
  gmi_fail("recompile with -DUSE_EGADS=ON");
}
#endif

struct gmi_model* gmi_egads_init(ego body, int nregions)
{
  eg_body = body;

  // set the context
  EG_getContext(eg_body, &eg_context);

  struct gmi_model *model;
  model = (struct gmi_model*)malloc(sizeof(*model));
  model->ops = &ops;

  int nverts, nedges, nfaces;

  EG_getBodyTopos(eg_body, NULL, NODE, &nverts, NULL);
  EG_getBodyTopos(eg_body, NULL, EDGE, &nedges, NULL);
  EG_getBodyTopos(eg_body, NULL, FACE, &nfaces, NULL);

  /// fix below to read supplementary file
  // I believe this should be shell, but always seems to result in 1 shell
  // EG_getBodyTopos(eg_body, NULL, SHELL, &nregions, NULL); // BODY?

  model->n[0] = nverts;
  model->n[1] = nedges;
  model->n[2] = nfaces;
  model->n[3] = nregions;

  for (int i = 0; i < 4; ++i)
  {
    egads_global_ents[i] = (egads_ent*)EG_alloc(sizeof(*egads_global_ents[i])
                                              * model->n[i]);
  }

  /// populate global array
  for (int dim = 0; dim < 4; ++dim)
  {
    for (int i = 0; i < model->n[dim]; ++i)
    {

      if (dim == 0)
      {
        EG_objectBodyTopo(eg_body, NODE, i+1,
                          &(egads_global_ents[dim][i].ego_ent));
      }
      else if (dim == 1)
      {
        EG_objectBodyTopo(eg_body, EDGE, i+1,
                          &(egads_global_ents[dim][i].ego_ent));
      }
      else if (dim == 2)
      {
        EG_objectBodyTopo(eg_body, FACE, i+1,
                          &(egads_global_ents[dim][i].ego_ent));
      }
      else if (dim == 3) // no EGADS 3D objects, just track with dim and tag
      {
        egads_global_ents[dim][i].ego_ent = NULL;
      }

      egads_global_ents[dim][i].dim = dim;
      egads_global_ents[dim][i].tag = i+1;
      if (dim < 3) {
        printf("ent dim: %d and tag: %d, magic number: %d\n", dim, i+1, (*egads_global_ents[dim][i].ego_ent).magicnumber);
        printf("actual dim: %d, tag: %d\n", egads_global_ents[dim][i].dim, egads_global_ents[dim][i].tag);
      }
    }
  }


  // printf("after:\n");
  // for (int dim = 0; dim < 4; ++dim)
  // {
  //   for (int i = 0; i < model->n[dim]; ++i)
  //   {
  //     if (dim < 3) {
  //       printf("ent dim: %d and tag: %d, magic number: %d\n", dim, i+1, (*egads_global_ents[dim][i].ego_ent).magicnumber);
  //       printf("actual dim: %d, tag: %d\n", egads_global_ents[dim][i].dim, egads_global_ents[dim][i].tag);
  //     }
  //   }
  // }

  return model;
}

void gmi_egads_init_adjacency(int ***adjacency)
{
  for (int i = 0; i < 6; ++i)
  {
    adjacency_table[i] = adjacency[i];
  }
}

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
