/****************************************************************************** 

  Copyright 2014 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/
#include <PCU.h>
#include "gmi_cap.h"
#include <stdlib.h>
#include <gmi.h>
#include <vector>
#include <pcu_util.h>


gmi_ent* toGmiEntity(M_GTopo topo)
{
  std::size_t hdl = topo.get();
  PCU_ALWAYS_ASSERT(hdl > 0);
  return reinterpret_cast<gmi_ent*>(hdl);
}

M_GTopo fromGmiEntity(gmi_ent* g)
{
  std::size_t hdl = reinterpret_cast<std::size_t>(g);
  M_GTopo topo;
  topo.set(hdl);
  return topo;
}


struct cap_model {
  struct gmi_model model;
  GeometryDatabaseInterface* geomInterface;
  bool owned;
};

static gmi_iter* begin(gmi_model* m, int dim)
{
  cap_model* cm = (cap_model*)m;
  M_GBRep brep;
  cm->geomInterface->get_brep_by_index(0, brep);
  GeometrySmartIterator* giter = new GeometrySmartIterator(cm->geomInterface);
  if (dim == 0)
    cm->geomInterface->get_topo_iterator(brep, VERTEX, *giter);
  if (dim == 1)
    cm->geomInterface->get_topo_iterator(brep, EDGE, *giter);
  if (dim == 2)
    cm->geomInterface->get_topo_iterator(brep, FACE, *giter);
  if (dim == 3)
    cm->geomInterface->get_topo_iterator(brep, REGION, *giter);
  cm->geomInterface->iterator_begin(*giter);
  return (gmi_iter*)(giter);
}

/* NOTE: giter is located at the first item in the list, therefore
 * gmi_next has to return it before calling iterator_next on giter
 */
static gmi_ent* next(gmi_model*m, gmi_iter* i)
{
  cap_model* cm = (cap_model*)m;
  GeometrySmartIterator* giter = (GeometrySmartIterator*)i;

  M_GTopo topo = cm->geomInterface->iterator_value(*giter);

  if (!cm->geomInterface->iterator_end(*giter))
    cm->geomInterface->iterator_next(*giter);
  else
    return 0;

  return toGmiEntity(topo);
}

static void end(gmi_model*, gmi_iter* i)
{
  GeometrySmartIterator* giter = (GeometrySmartIterator*)i;
  delete giter;
}

static int get_dim(gmi_model* m, gmi_ent* e)
{
  cap_model* cm = (cap_model*)m;
  M_GTopo topo = fromGmiEntity(e);
  if (cm->geomInterface->is_vertex(topo))
    return 0;
  if (cm->geomInterface->is_edge(topo))
    return 1;
  if (cm->geomInterface->is_face(topo))
    return 2;
  if (cm->geomInterface->is_region(topo))
    return 3;
  return 0;
}

static int get_tag(gmi_model* m, gmi_ent* e)
{
  cap_model* cm = (cap_model*)m;
  M_GTopo topo = fromGmiEntity(e);
  std::size_t id;
  cm->geomInterface->get_id(topo, id);
  return (int)id;
}

static gmi_ent* find(gmi_model* m, int dim, int tag)
{
  (void)m;
  (void)dim;
  (void)tag;
  printf("_find_ not implemented!\n");
  return 0;
}

static gmi_set* adjacent(gmi_model* m, gmi_ent* e, int dim)
{
  (void)m;
  (void)e;
  (void)dim;
  printf("_adjacent_ not implemented!\n");
  return 0;
}

static void eval(struct gmi_model* m, struct gmi_ent* e,
      double const p[2], double x[3])
{
  cap_model* cm = (cap_model*)m;
  M_GTopo topo = fromGmiEntity(e);
  vec3d point;
  cm->geomInterface->get_point(topo, vec3d(p[0], p[1], 0.0), point);
  x[0] = point[0];
  x[1] = point[1];
  x[2] = point[2];
}

static void reparam(struct gmi_model* m, struct gmi_ent* from,
      double const from_p[2], struct gmi_ent* to, double to_p[2])
{
  cap_model* cm = (cap_model*)m;
  M_GTopo fromTopo = fromGmiEntity(from);
  M_GTopo toTopo   = fromGmiEntity(to);
  vec3d fromParam(from_p[0], from_p[1], 0.0);
  vec3d toParam;
  cm->geomInterface->reparametrize(fromTopo, fromParam, toTopo, toParam);
  to_p[0] = toParam[0];
  to_p[1] = toParam[1];
}

static int periodic(struct gmi_model* m, struct gmi_ent* e, int dim)
{
  cap_model* cm = (cap_model*)m;
  M_GTopo topo = fromGmiEntity(e);
  int paramType;
  cm->geomInterface->get_parametrization_info(topo, dim, paramType);
  // paramType is a bit mask with the following definitions
  // the bit value at location 1 determines periodicity
  /* enum ParametrizationType */
  /* { */
  /*     PARAM_UNKNOWN=0,      //!< Unknown parametrization */
  /*     PARAM_CONTINUOUS=1,   //!< Continuous parametrization */
  /*     PARAM_PERIODIC=2,     //!< Periodic parametrization */
  /*     PARAM_COMPOSITE=4,    //!< parametrization of a composite edge or face */
  /*     PARAM_UNBOUNDED=8     //!< infinite parametrization */
  /* }; */
  return paramType & (1<<1);
}

static void range(struct gmi_model* m, struct gmi_ent* e, int dim,
    double r[2])
{
  cap_model* cm = (cap_model*)m;
  M_GTopo topo = fromGmiEntity(e);
  double lower, upper;
  cm->geomInterface->get_parametrization_range(topo, dim, lower, upper);
  r[0] = lower;
  r[1] = upper;
}

static void closest_point(struct gmi_model* m, struct gmi_ent* e,
    double const from[3], double to[3], double to_p[2])
{
  (void)m;
  (void)from;
  (void)e;
  (void)to;
  (void)to_p;
  printf("_closest_point_ not implemented!\n");
}

static void normal(struct gmi_model* m, struct gmi_ent* e,
    double const p[2], double n[3])
{
  (void)m;
  (void)e;
  (void)p;
  (void)n;
  printf("_normal_ not implemented!\n");
}

static void first_derivative(struct gmi_model* m, struct gmi_ent* e,
    double const p[2], double t0[3], double t1[3])
{
  cap_model* cm = (cap_model*)m;
  M_GTopo topo = fromGmiEntity(e);
  vec3d param(p[0], p[1], 0.0);
  std::vector<double> dxyz;
  cm->geomInterface->get_derivative(topo, param, 1, dxyz);
  t0[0] = dxyz[0];
  t0[1] = dxyz[1];
  t0[2] = dxyz[2];
  t1[0] = dxyz[3];
  t1[1] = dxyz[4];
  t1[2] = dxyz[5];
}

static int is_point_in_region(struct gmi_model* m, struct gmi_ent* e,
    double point[3])
{
  (void)m;
  (void)e;
  (void)point;
  printf("_is_point_in_region_ not implemented!\n");
  return 0;
}

static int is_in_closure_of(struct gmi_model* m, struct gmi_ent* e,
    struct gmi_ent* et)
{
  (void)m;
  (void)e;
  (void)et;
  printf("_is_in_closure_of_ not implemented!\n");
  return 0;
}

static int is_discrete_ent(struct gmi_model*, struct gmi_ent* e)
{
  (void)e;
  printf("_is_discrete_ent_ not implemented!\n");
  return 0;
}

static void destroy(gmi_model* m)
{
  (void)m;
  printf("_destroy_ not implemented!\n");
}

static struct gmi_model_ops ops;


void gmi_cap_start(void)
{
}

void gmi_cap_stop(void)
{
}

void gmi_register_cap(void)
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
  ops.is_discrete_ent = is_discrete_ent;
  ops.destroy = destroy;
  /* gmi_register(create_smd, "smd"); */
  /* gmi_register(create_native, "xmt_txt"); */
  /* gmi_register(create_native, "x_t"); */
  /* gmi_register(create_native, "sat"); */
}

/* static gmi_model* owned_import(GeometryDatabaseInterface* gi) */
/* { */
/*   (void)gi; */
/*   printf("not implemented!\n"); */
/*   return 0; */
/* } */

gmi_model* gmi_import_cap(GeometryDatabaseInterface* gi)
{
  cap_model* m;
  m = (cap_model*)malloc(sizeof(*m));
  m->model.ops = &ops;
  m->geomInterface = gi;
  M_GBRep brep;
  int numBreps;

  m->geomInterface->get_num_breps(numBreps);
  PCU_ALWAYS_ASSERT(numBreps == 1);
  m->geomInterface->get_brep_by_index(0, brep);
  m->geomInterface->get_num_topos(brep, VERTEX, m->model.n[0]);
  m->geomInterface->get_num_topos(brep, EDGE,	m->model.n[1]);
  m->geomInterface->get_num_topos(brep, FACE,	m->model.n[2]);
  m->geomInterface->get_num_topos(brep, REGION, m->model.n[3]);
  m->owned = false;
  return &m->model;
}

GeometryDatabaseInterface* gmi_export_cap(gmi_model* m)
{
  cap_model* cm = (cap_model*)m;
  return cm->geomInterface;
}
