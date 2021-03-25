/******************************************************************************

  Copyright 2014 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/
#include <dlfcn.h>
#include <math.h>
#include <string.h>
#include <stdio.h>

#include "egads.h"

#include "PCU.h"
#include "gmi.h"
#include "gmi_egads.h"
#include "gmi_egads_config.h"

/**
 * The following are the EGADS APIs we use to implement the GMI interface
 * Since we'll load EGADS or EGADSlite dynamically at runtime each of these
 * functions needs to be a function pointer loaded with dlsym
*/

/// handle to dynamically loaded EGADS/EGADSlite library
static void *handle;

/// EGADS API of dynamically loaded functions
/* memory functions */
static void * (*EGADS_alloc)(size_t nbytes);
static void (*EGADS_free)(void *pointer);
/* base-level object functions */
static int (*EGADS_open)( ego *context );
static int (*EGADS_loadModel)(ego context, int bflg, const char *name, 
                              ego *model);
static int (*EGADS_getContext)(ego object, ego *context);
static int (*EGADS_close)(ego context);
/* geometry functions */
static int (*EGADS_getRange)(const ego geom, double *range, int *periodic);
static int (*EGADS_evaluate)(const ego geom, /*@null@*/ const double *param, 
                             double *results);
static int (*EGADS_invEvaluate)(const ego geom, double *xyz, double *param,
                                double *results);
/* topology functions */
static int (*EGADS_getTopology)(const ego topo, ego *geom, int *oclass, 
                                int *type, /*@null@*/ double *limits, 
                                int *nChildren, ego **children, int **sense);
static int (*EGADS_getBodyTopos)(const ego body, /*@null@*/ ego src,
                                 int oclass, int *ntopo,
                                 /*@null@*/ ego **topos);
static int (*EGADS_indexBodyTopo)(const ego body, const ego src);
static int (*EGADS_objectBodyTopo)(const ego body, int oclass, int index,
                                   ego *obj);
static int (*EGADS_inTopology)(const ego topo, const double *xyz);
static int (*EGADS_getEdgeUV)(const ego face, const ego edge, int sense,
                              double t, double *UV);
static int (*EGADS_getBoundingBox)(const ego topo, double *bbox);
static int (*EGADS_isEquivalent)(const ego topo1, const ego topo2);

/** \brief initialize a gmi_model with an EGADS body and number of regions */
struct gmi_model* gmi_egads_init(struct egObject *body, int numRegions);

/** \brief initialize the model adjacency table for 3D regions */
void gmi_egads_init_adjacency(struct gmi_model* m, int ***adjacency);

/// global context used by EGADS initialized by `gmi_egads_start`
static ego eg_context;

typedef struct egads_ent
{
  ego ego_ent;
  int dim;
  int tag;
} egads_ent;

typedef struct egads_model
{
  struct gmi_model model;
  ego eg_body;

  /// array of model entities
  /// egads_model_ents[0] - array of model verts
  /// egads_model_ents[1] - array of model edges
  /// egads_model_ents[2] - array of model faces
  /// egads_model_ents[3] - array of model regions
  egads_ent *egads_model_ents[4];

  /// adjacency table to be used for the "dummy" 3D elements that EGADS doesn't
  /// natively track
  int **adjacency_table[6];
} egads_model;

typedef struct egads_iter
{
  egads_ent *ents;
  int nelem;
  int idx;
} egads_iter;

/// read adjacency table from binary file
static void read_adj_table(const char* filename,
                           int *nregions,
                           int *** adjacency_table)
{
  /// filename of supplementary model file
  /// 4 chars longer for ".sup" suffix plus 1 for string termination character
  char *sup_filename;
  sup_filename = EGADS_alloc(strlen(filename)+4+1);
  if (sup_filename == NULL)
  {
    gmi_fail("failed to allocate memory for new string");
  }
  sup_filename[0] = '\0';
  strcat(sup_filename, filename);
  strcat(sup_filename, ".sup");

  FILE *adj_table_file = fopen(sup_filename, "rb");
  EGADS_free(sup_filename);
  if (adj_table_file == NULL)
  {
    gmi_fail("failed to open supplementary EGADS model file!");
  }

  int header[6];
  size_t count = fread(header, sizeof(int), 6, adj_table_file);
  if (count != 6) gmi_fail("fread failed!\n");

  *nregions = header[0];

  for (int i = 0; i < 6; ++i)
  {
    adjacency_table[i] = (int**)EGADS_alloc(sizeof(*(adjacency_table[i]))
                                         * header[i]);
    if (adjacency_table[i] == NULL) {
      char fail[50];
      sprintf(fail, "failed to alloc memory for adjacency_table[%d]", i);
      /// TODO: this could cause memory leak
      gmi_fail(fail);
    }
    for (int j = 0; j < header[i]; ++j) {
      int nadjacent = -1;
      count = fread(&nadjacent, sizeof(int), 1, adj_table_file);
      if (count != 1) gmi_fail("fread failed!\n"); 
      adjacency_table[i][j] = (int*)EGADS_alloc(sizeof(*(adjacency_table[i][j]))
                                             * (nadjacent+1));
      if (adjacency_table[i][j] == NULL) {
        char fail[100];
        sprintf(fail, "failed to alloc memory for "
                "adjacency_table[%d][%d]", i,j);
        /// TODO: this could cause memory leak
        gmi_fail(fail);
      }
      adjacency_table[i][j][0] = nadjacent;
      count = fread(&(adjacency_table[i][j][1]), sizeof(int), nadjacent,
                    adj_table_file);
      if (count != (size_t)nadjacent) gmi_fail("fread failed!\n");
    }
  }
  fclose(adj_table_file);
}


/// TODO: consider optimizing adjacency tables and access
static void get_3D_adjacency(struct gmi_model* m,
                             egads_ent* ent, 
                             int adj_dim, 
                             int *num_adjacent, 
                             egads_ent*** adjacent_ents)
{
  egads_model *egm = (egads_model*)m;

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

  int *adj_tags = egm->adjacency_table[pairing][ent_tag-1];
  *num_adjacent = adj_tags[0];
  *adjacent_ents = (egads_ent**)EGADS_alloc(sizeof(**adjacent_ents) * (*num_adjacent));

  for (int i = 0; i < *num_adjacent; ++i)
  {
    (*adjacent_ents)[i] = (egads_ent*)m->ops->find(m, adj_dim, adj_tags[i+1]);
  }
}

/// TODO: implement based on adjacent face's bounding boxes?
static void get_3D_bounding_box(egads_ent *ent, double *box)
{
  (void)ent;
  (void)box;
  gmi_fail("3D bounding box not implemented!\n");
}

/// reparameterize a vertex onto an edge
static void getVertexT(struct gmi_model* m,
                       struct gmi_ent* to,
                       struct gmi_ent* from,
                       double* t)
{
  double diff;
  double t_range[2];
  m->ops->range(m, to, 0, &(t_range[0]));
  double vtx_pnt[3];
  double p[] = {0, 0};
  m->ops->eval(m, from, p, &(vtx_pnt[0]));
  double t_pnt[3];
  m->ops->eval(m, to, &(t_range[0]), &(t_pnt[0]));
  diff = sqrt(pow(vtx_pnt[0] - t_pnt[0], 2) + 
              pow(vtx_pnt[1] - t_pnt[1], 2) + 
              pow(vtx_pnt[2] - t_pnt[2], 2));
  if (diff < 0.001)
  {
    *t = t_range[0];
  }
  else
  {
    m->ops->eval(m, to, &(t_range[1]), &(t_pnt[0]));
    // diff = sqrt(pow(vtx_pnt[0] - t_pnt[0], 2) + 
    //             pow(vtx_pnt[1] - t_pnt[1], 2) + 
    //             pow(vtx_pnt[2] - t_pnt[2], 2));
    *t = t_range[1];
  }
  return;
}

/// TODO: handle degenerate edges (cant eval on them)
/// function to reparameterize a vertex onto a face
/// (get the u,v parametric coodinates associated with the vertex)
static void getVertexUV(struct gmi_model* m,
                        struct gmi_ent* to,
                        struct gmi_ent* from,
                        double to_p[2])
{
  struct gmi_set* adj_faces = NULL;
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
    if (adj_faces != NULL)
      gmi_free_set(adj_faces);
    gmi_free_set(adj_edges);
    return;
}

static struct gmi_iter* begin(struct gmi_model* m, int dim)
{
  egads_model *egm = (egads_model*)m;

  struct egads_iter *eg_iter;
  if (dim >= 0 && dim <= 3)
  {
    eg_iter = EGADS_alloc(sizeof(*eg_iter));
    if (eg_iter == NULL)
    {
      gmi_fail("EGADS_alloc failed to allocate memory for iter");
      return (struct gmi_iter*)NULL;
    }
    int nents = egm->model.n[dim];
    eg_iter->ents = &(egm->egads_model_ents[dim][0]);
    eg_iter->nelem = nents;
    eg_iter->idx = 0;
    return (struct gmi_iter*)eg_iter;
  }
  return (struct gmi_iter*)NULL;
}

static struct gmi_ent* next(struct gmi_model* m, struct gmi_iter* i)
{
  (void)m;
  struct egads_iter *eg_iter = (struct egads_iter*)i;
  if (eg_iter->idx < eg_iter->nelem)
  {
    return (struct gmi_ent*)(eg_iter->ents+eg_iter->idx++);
  }
  return (struct gmi_ent*)NULL;
}

static void end(struct gmi_model* m, struct gmi_iter* i)
{
  (void)m;

  struct egads_iter *eg_iter = (struct egads_iter*)i;
  /// I think freeing the array here will free it too early, 
  if (eg_iter != NULL)
  {
    EGADS_free(eg_iter);
  }
}

static int get_dim(struct gmi_model* m, struct gmi_ent* e)
{
  (void)m;
  egads_ent* eg_ent = (egads_ent*)e;
  return eg_ent->dim;
}

static int get_tag(struct gmi_model* m, struct gmi_ent* e)
{
  (void)m;
  egads_ent* eg_ent = (egads_ent*)e;
  return eg_ent->tag;
}

static struct gmi_ent* find(struct gmi_model* m, int dim, int tag)
{
  egads_model *egm = (egads_model*)m;
  return (struct gmi_ent*)&(egm->egads_model_ents[dim][tag-1]);
}

// find (dim)-dimensional entities adjacent to e
static struct gmi_set* adjacent(struct gmi_model* m, 
                                struct gmi_ent* e, 
                                int dim)
{
  egads_model *egm = (egads_model*)m;

  egads_ent *eg_ent = (egads_ent*)e;

  int num_adjacent = 0;

  egads_ent **adjacent_ents = NULL;

  if (eg_ent->dim == 3 || dim == 3)
  {
    get_3D_adjacency(m, eg_ent, dim, &num_adjacent, &adjacent_ents);
  }
  else // only dealing with egos
  {
    ego *adjacent_egos;
    if (dim == 0)
      EGADS_getBodyTopos(egm->eg_body, (eg_ent->ego_ent), 20, &num_adjacent, &adjacent_egos);
    else if (dim == 1)
      EGADS_getBodyTopos(egm->eg_body, (eg_ent->ego_ent), 21, &num_adjacent, &adjacent_egos);
    else if (dim == 2)
      EGADS_getBodyTopos(egm->eg_body, (eg_ent->ego_ent), 23, &num_adjacent, &adjacent_egos);

    adjacent_ents = (egads_ent**)EGADS_alloc(sizeof(*adjacent_ents) * num_adjacent);
    for (int i = 0; i < num_adjacent; i++)
    {
      int adj_ent_tag = EGADS_indexBodyTopo(egm->eg_body, adjacent_egos[i]);
      adjacent_ents[i] = (egads_ent*)m->ops->find(m, dim, adj_ent_tag);
    }
    EGADS_free(adjacent_egos);
  }

  struct gmi_set *gmi_adj_ent = gmi_make_set(num_adjacent);
  for (int i = 0; i < num_adjacent; ++i)
  {
    gmi_adj_ent->e[i] = (struct gmi_ent*)adjacent_ents[i];
  }
  EGADS_free(adjacent_ents);
  adjacent_ents = NULL;

  return gmi_adj_ent;
}

/// TODO: error check to make sure ego_ent != NULL?
static void eval(struct gmi_model* m, 
                 struct gmi_ent* e,
                 double const p[2],
                 double x[3])
{
  double results[18];
  egads_ent *eg_ent = (egads_ent*)e;
  ego ego_ent = eg_ent->ego_ent;
  int dim = m->ops->dim(m, e);
  if (dim > 0 && dim < 3)
  {
    EGADS_evaluate(ego_ent, p, results);
    x[0] = results[0];
    x[1] = results[1];
    x[2] = results[2];
  }
  else if (dim == 0)
  {
    double data[4];
    int oclass, mtype, nbody, *senses;
    ego geom, *eg_bodies;
    EGADS_getTopology(ego_ent, &geom, &oclass, &mtype, data, &nbody,
                   &eg_bodies, &senses);
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
static void reparam(struct gmi_model* m,
                    struct gmi_ent* from,
                    double const from_p[2],
                    struct gmi_ent* to,
                    double to_p[2])
{
  int from_dim, to_dim;
  from_dim = get_dim(m, from);
  to_dim = get_dim(m, to);
  egads_ent *eg_ent_from = (egads_ent*)from;
  ego ego_from = eg_ent_from->ego_ent;

  egads_ent *eg_ent_to = (egads_ent*)to;
  ego ego_to = eg_ent_to->ego_ent;

  if ((from_dim == 1) && (to_dim == 2))
  {
    EGADS_getEdgeUV(ego_to, ego_from, 1, from_p[0], to_p);
    return;
  }
  if ((from_dim == 0) && (to_dim == 2))
  {
    getVertexUV(m, to, from, to_p);
    return;
  }
  if ((from_dim == 0) && (to_dim == 1))
  {
    getVertexT(m, to, from, &to_p[0]);
    return;
  }
  gmi_fail("bad dimensions in gmi_egads reparam");
}

/// TODO: error check to make sure ego_ent != NULL?
static int periodic(struct gmi_model* m,
                    struct gmi_ent* e,
                    int dir)
{
  int ent_dim = get_dim(m, e);
  int periodic;
  egads_ent *eg_ent = (egads_ent*)e;
  ego ego_ent = eg_ent->ego_ent;

  double range[4];
  EGADS_getRange(ego_ent, range, &periodic);

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
static void range(struct gmi_model* m,
                  struct gmi_ent* e,
                  int dir,
                  double r[2])
{
  int ent_dim = get_dim(m, e);
  double range[4];
  int periodic;
  egads_ent *eg_ent = (egads_ent*)e;
  ego ego_ent = eg_ent->ego_ent;

  EGADS_getRange(ego_ent, range, &periodic);
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
static void closest_point(struct gmi_model* m,
                          struct gmi_ent* e, 
                          double const from[3],
                          double to[3],
                          double to_p[2])
{
  (void)m;
  egads_ent *eg_ent = (egads_ent*)e;
  ego ego_ent = eg_ent->ego_ent;
  double xyz[] = {from[0], from[1], from[2]};
  EGADS_invEvaluate(ego_ent, &xyz[0], &to_p[0], &to[0]);
}

/// TODO: error check to make sure ego_ent != NULL?
static void normal(struct gmi_model* m,
                   struct gmi_ent* e,
                   double const p[2],
                   double n[3])
{
  double du[3], dv[3];
  m->ops->first_derivative(m, e, p, du, dv);
  // cross du and dv to get n
  n[0] = du[1]*dv[2] - du[2]*dv[1];
  n[1] = du[2]*dv[0] - du[0]*dv[2];
  n[2] = du[0]*dv[1] - du[1]*dv[0];

  egads_ent *eg_ent = (egads_ent*)e;
  ego ego_ent = eg_ent->ego_ent;
  double data[4];
  int oclass, mtype, nbody, *senses;
  ego geom, *eg_bodies;
  EGADS_getTopology(ego_ent, &geom, &oclass, &mtype, data, &nbody, &eg_bodies, &senses);

  double n_mag = sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
  n[0] *= mtype / n_mag;
  n[1] *= mtype / n_mag;
  n[2] *= mtype / n_mag;
}

/// TODO: error check to make sure ego_ent != NULL?
static void first_derivative(struct gmi_model* m,
                             struct gmi_ent* e,
                             double const p[2],
                             double t0[3],
                             double t1[3])
{
  int ent_dim = get_dim(m, e);
  double results[18];
  egads_ent *eg_ent = (egads_ent*)e;
  ego ego_ent = eg_ent->ego_ent;
  EGADS_evaluate(ego_ent, p, results);
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

/// TODO: make this work for 3D entity
static int is_point_in_region(struct gmi_model* m,
                              struct gmi_ent* e,
                              double p[3])
{
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
    int status = EGADS_inTopology(ego_ent, p);
    if (status == EGADS_SUCCESS)
      return 1;
    else
      return 0;
  }
}

/// TODO: make this work for new 3D object
static void bbox(struct gmi_model* m,
                 struct gmi_ent* e,
                 double bmin[3],
                 double bmax[3])
{
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
    EGADS_getBoundingBox(ego_ent, box);
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
static int is_in_closure_of(struct gmi_model* m,
                            struct gmi_ent* e,
                            struct gmi_ent* et)
{
  egads_model *egm = (egads_model*)m;

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
        EGADS_free(adjacent_ents);
        return 1;
      }
    }
    EGADS_free(adjacent_ents);
    return 0;
  }
  else
  {
    ego *adjacent_egos = NULL;
    if (ent_dim == 0)
      EGADS_getBodyTopos(egm->eg_body, ego_region, NODE, &num_adjacent, &adjacent_egos);
    else if (ent_dim == 1)
      EGADS_getBodyTopos(egm->eg_body, ego_region, EDGE, &num_adjacent, &adjacent_egos);
    else if (ent_dim == 2)
      EGADS_getBodyTopos(egm->eg_body, ego_region, FACE, &num_adjacent, &adjacent_egos);
    for (int i = 0; i < num_adjacent; ++i)
    {
      if (EGADS_isEquivalent(ego_ent, adjacent_egos[i]))
      {
        EGADS_free(adjacent_egos);
        return 1;
      }
    }
    EGADS_free(adjacent_egos);
    return 0;
  }
}

/// what does this function do?
static int is_discrete_ent(struct gmi_model* m, struct gmi_ent* e)
{
  (void)m;
  (void)e;
  gmi_fail("is_discrete_ent not implemented");
}

static void destroy(struct gmi_model* m)
{
  egads_model *egm = (egads_model*)m;

  for (int dim = 0; dim < 4; ++dim)
  {
    EGADS_free(egm->egads_model_ents[dim]);
  }

  int sizes[] = {m->n[3], m->n[3], m->n[3],
                 m->n[0], m->n[1], m->n[2]};
  for (int i = 0; i < 6; ++i)
  {
    for (int j = 0; j < sizes[i]; ++j)
    {
      EGADS_free(egm->adjacency_table[i][j]);
    }
    EGADS_free(egm->adjacency_table[i]);
  }

  free(egm);
}

static struct gmi_model_ops ops;

/// TODO: Come up with a better flag? 
/// TODO: re-write for EGADSlite - model loading is different
// #ifdef HAVE_EGADS
#if 1
struct gmi_model* gmi_egads_load(const char* filename)
{
  ego eg_model;
  int load_status = EGADS_loadModel(eg_context, 0, filename, &eg_model);
  if (load_status != EGADS_SUCCESS)
  {
    char str[50]; // big enough
    sprintf(str, "EGADS failed to load model with error code: %d", load_status);
    gmi_fail(str);
  }

  /// TODO: only store the outputs I need, replace the rest with NULL
  int oclass, mtype, nbody, *senses;
  ego geom, *eg_bodies;
  int status = EGADS_getTopology(eg_model, &geom, &oclass, &mtype, NULL, &nbody,
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
  int **adjacency_table[6];
  read_adj_table(filename, &nregions, adjacency_table);

  struct gmi_model* m = gmi_egads_init(eg_bodies[0], nregions);
  gmi_egads_init_adjacency(m, adjacency_table);
  return m;
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
  // set the context
  EGADS_getContext(body, &eg_context);

  egads_model* m;
  m = (egads_model*)EGADS_alloc(sizeof(*m));
  m->model.ops = &ops;
  m->eg_body = body;

  int nverts, nedges, nfaces;
  EGADS_getBodyTopos(body, NULL, NODE, &nverts, NULL);
  EGADS_getBodyTopos(body, NULL, EDGE, &nedges, NULL);
  EGADS_getBodyTopos(body, NULL, FACE, &nfaces, NULL);

  m->model.n[0] = nverts;
  m->model.n[1] = nedges;
  m->model.n[2] = nfaces;
  m->model.n[3] = nregions;

  for (int i = 0; i < 4; ++i)
  {
    m->egads_model_ents[i] = (egads_ent*)EGADS_alloc(sizeof(*m->egads_model_ents[i])
                                              * m->model.n[i]);
  }

  /// populate model entity array
  for (int dim = 0; dim < 4; ++dim)
  {
    for (int i = 0; i < m->model.n[dim]; ++i)
    {
      if (dim == 0)
      {
        EGADS_objectBodyTopo(body, NODE, i+1,
                          &(m->egads_model_ents[dim][i].ego_ent));
      }
      else if (dim == 1)
      {
        EGADS_objectBodyTopo(body, EDGE, i+1,
                          &(m->egads_model_ents[dim][i].ego_ent));
      }
      else if (dim == 2)
      {
        EGADS_objectBodyTopo(body, FACE, i+1,
                          &(m->egads_model_ents[dim][i].ego_ent));
      }
      else if (dim == 3) // no EGADS 3D objects, just track with dim and tag
      {
        m->egads_model_ents[dim][i].ego_ent = NULL;
      }

      m->egads_model_ents[dim][i].dim = dim;
      m->egads_model_ents[dim][i].tag = i+1;
    }
  }

  return &m->model;
}

void gmi_egads_init_adjacency(struct gmi_model* m, int ***adjacency)
{
  egads_model *egm = (egads_model*)m;
  for (int i = 0; i < 6; ++i)
  {
    egm->adjacency_table[i] = adjacency[i];
  }
}

static void open_egads_lib(void)
{
  const int rank = PCU_Proc_Self();
  if (rank == 0)
  {
    handle = dlopen(EGADS_LIBRARIES, RTLD_NOW);
  }
  else
  {
    handle = dlopen(EGADSLITE_LIBRARIES, RTLD_NOW);
  }
}

static void close_egads_lib(void)
{
  dlclose(handle);
}

static void load_dl_function(void (**fptr)(void), const char *symbol)
{
  *fptr = (void (*)(void))dlsym(handle, symbol);
  if (!*fptr) {
    /* no such symbol */
    printf("Symbol %s not found!\n", symbol);
    fprintf(stderr, "Error: %s\n", dlerror());
    dlclose(handle);
    char str[100]; // big enough
    sprintf(str, "Symbol %s not found!\n", symbol);
    gmi_fail(str);
  }
}

static void init_egads_functions(void)
{
  load_dl_function((void (**)(void))&EGADS_alloc, "EG_alloc");
  load_dl_function((void (**)(void))&EGADS_free, "EG_free");
  load_dl_function((void (**)(void))&EGADS_open, "EG_open");
  load_dl_function((void (**)(void))&EGADS_loadModel, "EG_loadModel");
  load_dl_function((void (**)(void))&EGADS_getContext, "EG_getContext");
  load_dl_function((void (**)(void))&EGADS_close, "EG_close");
  load_dl_function((void (**)(void))&EGADS_getRange, "EG_getRange");
  load_dl_function((void (**)(void))&EGADS_evaluate, "EG_evaluate");
  load_dl_function((void (**)(void))&EGADS_invEvaluate, "EG_invEvaluate");
  load_dl_function((void (**)(void))&EGADS_getTopology, "EG_getTopology");
  load_dl_function((void (**)(void))&EGADS_getBodyTopos, "EG_getBodyTopos");
  load_dl_function((void (**)(void))&EGADS_indexBodyTopo, "EG_indexBodyTopo");
  load_dl_function((void (**)(void))&EGADS_objectBodyTopo, "EG_objectBodyTopo");
  load_dl_function((void (**)(void))&EGADS_inTopology, "EG_inTopology");
  load_dl_function((void (**)(void))&EGADS_getEdgeUV, "EG_getEdgeUV");
  load_dl_function((void (**)(void))&EGADS_getBoundingBox, "EG_getBoundingBox");
  load_dl_function((void (**)(void))&EGADS_isEquivalent, "EG_isEquivalent");

  // *(void **)(&EGADS_alloc) = dlsym(handle, "EG_alloc");
  // *(void **)(&EGADS_free) = dlsym(handle, "EG_free");
  // *(void **)(&EGADS_open) = dlsym(handle, "EG_open");
  // *(void **)(&EGADS_loadModel) = dlsym(handle, "EG_loadModel");
  // *(void **)(&EGADS_getContext) = dlsym(handle, "EG_getContext");
  // *(void **)(&EGADS_close) = dlsym(handle, "EG_close");
  // *(void **)(&EGADS_getRange) = dlsym(handle, "EG_getRange");
  // *(void **)(&EGADS_evaluate) = dlsym(handle, "EG_evaluate");
  // *(void **)(&EGADS_invEvaluate) = dlsym(handle, "EG_invEvaluate");
  // *(void **)(&EGADS_getTopology) = dlsym(handle, "EG_getTopology");
  // *(void **)(&EGADS_getBodyTopos) = dlsym(handle, "EG_getBodyTopos");
  // *(void **)(&EGADS_indexBodyTopo) = dlsym(handle, "EG_indexBodyTopo");
  // *(void **)(&EGADS_objectBodyTopo) = dlsym(handle, "EG_objectBodyTopo");
  // *(void **)(&EGADS_inTopology) = dlsym(handle, "EG_inTopology");
  // *(void **)(&EGADS_getEdgeUV) = dlsym(handle, "EG_getEdgeUV");
  // *(void **)(&EGADS_getBoundingBox) = dlsym(handle, "EG_getBoundingBox");
  // *(void **)(&EGADS_isEquivalent) = dlsym(handle, "EG_isEquivalent");
}

void gmi_egads_start(void)
{
  
  open_egads_lib();
  init_egads_functions();

  int status = EGADS_open(&eg_context);
  if (status != EGADS_SUCCESS)
  {
    char str[50]; // big enough
    sprintf(str, "EGADS failed to open with error code: %d", status);
    gmi_fail(str);
  }
}

void gmi_egads_stop(void)
{
  EGADS_close(eg_context);
  close_egads_lib();
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
