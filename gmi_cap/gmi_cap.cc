/****************************************************************************** 

  Copyright 2014 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/
#include "gmi_cap.h"
#include <stdlib.h>
#include <gmi.h>
#include <vector>
#include <pcu_util.h>
#include <lionPrint.h>


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
  gmi_fail("unable to determine model dimension\n");
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
  cap_model* cm = (cap_model*)m;
  M_GTopo topo;
  if (dim == 0)
    topo = cm->geomInterface->get_topo_by_id(VERTEX, tag);
  else if (dim == 1)
    topo = cm->geomInterface->get_topo_by_id(EDGE, tag);
  else if (dim == 2)
    topo = cm->geomInterface->get_topo_by_id(FACE, tag);
  else if (dim == 3)
    topo = cm->geomInterface->get_topo_by_id(REGION, tag);
  else
    gmi_fail("input dim is out of range.");
  return toGmiEntity(topo);
}

static gmi_set* vector_to_set(std::vector<M_GTopo> gtopos)
{
  gmi_set* s = gmi_make_set(gtopos.size());
  for (int i = 0; i < s->n; i++) {
    s->e[i] = toGmiEntity(gtopos[i]);
  }
  return s;
}

static gmi_set* adjacent(gmi_model* m, gmi_ent* e, int dim)
{
  cap_model* cm = (cap_model*)m;
  M_GTopo gtopo = fromGmiEntity(e);
  int edim = gmi_dim(m, e);
  std::vector<M_GTopo> gtopos;
  if (edim == 0 && dim == 1)
    cm->geomInterface->get_adjacency(gtopo, EDGE, gtopos);
  else if (edim == 1 && dim == 0)
    cm->geomInterface->get_adjacency(gtopo, VERTEX, gtopos);
  else if (edim == 1 && dim == 2)
    cm->geomInterface->get_adjacency(gtopo, FACE, gtopos);
  else if (edim == 2 && dim == 0)
    cm->geomInterface->get_adjacency(gtopo, VERTEX, gtopos);
  else if (edim == 2 && dim == 1)
    cm->geomInterface->get_adjacency(gtopo, EDGE, gtopos);
  else if (edim == 2 && dim == 3)
    cm->geomInterface->get_adjacency(gtopo, REGION, gtopos);
  else if (edim == 1 && dim == 3)
    cm->geomInterface->get_adjacency(gtopo, REGION, gtopos);
  else if (edim == 3 && dim == 2)
    cm->geomInterface->get_adjacency(gtopo, FACE, gtopos);
  else
    gmi_fail("requested adjacency not available\n");
  return vector_to_set(gtopos);
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
  cap_model* cm = reinterpret_cast<cap_model*>(m);
  M_GTopo topo = fromGmiEntity(e);
  vec3d xyz_in, xyz_out, param_out;
  xyz_in[0] = from[0]; xyz_in[1] = from[1]; xyz_in[2] = from[2];
  // Capstone recommends replacing nullptr with `seedparam` -- a known nearby
  // parametric point. This interface doesn't have that functionality, but if
  // it becomes added, this call should be updated.
  MG_API_CALL(cm->geomInterface, get_closest_point_param(topo, xyz_in, nullptr,
    xyz_out, param_out));
  to[0] = xyz_out[0]; to[1] = xyz_out[1]; to[2] = xyz_out[2];
  to_p[0] = param_out[0]; to_p[1] = param_out[1];
}

static void normal(struct gmi_model* m, struct gmi_ent* e,
    double const p[2], double n[3])
{
  cap_model* cm = (cap_model*)m;
  M_GTopo topo = fromGmiEntity(e);
  vec2d param(p[0], p[1]);
  vec3d norm;
  cm->geomInterface->get_face_normal_parametrization(topo, param, norm);
  n[0] = norm[0];
  n[1] = norm[1];
  n[2] = norm[2];
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
  cap_model* cm = (cap_model*) m;
  std::vector<M_GTopo> topos;
  MG_API_CALL(cm->geomInterface, find_point_containment(vec3d(point),
    Geometry::REGION, topos, 0.0));
  M_GTopo gtopo = fromGmiEntity(e);
  for (size_t i = 0; i < topos.size(); ++i) {
    if (topos[i] == gtopo) return 1;
  }
  return 0;
}

static int is_in_closure_of(struct gmi_model* m, struct gmi_ent* e,
    struct gmi_ent* et)
{
  (void)m;
  (void)e;
  (void)et;
  if (get_dim(m, e) < get_dim(m, et)) {
    cap_model* cm = (cap_model*)m;
    M_GTopo ce = fromGmiEntity(e);
    M_GTopo cet = fromGmiEntity(et);
    GeometryTopo ce_type;
    MG_API_CALL(cm->geomInterface, get_topo_type(ce, ce_type));
    std::vector<M_GTopo> adj;
    MG_API_CALL(cm->geomInterface, get_adjacency(cet, ce_type, adj));
    for (size_t i = 0; i < adj.size(); ++i) {
      if (ce == adj[i]) {
        return 1;
      }
    }
  }
  return 0;
}

static void bbox(struct gmi_model* m, struct gmi_ent* e,
    double bmin[3], double bmax[3])
{
  cap_model* cm = (cap_model*)m;
  M_GTopo topo = fromGmiEntity(e);
  BBox bbox;
  cm->geomInterface->get_bounding_box(topo, bbox);
  vec3d bmi = bbox.min_value();
  vec3d bma = bbox.max_value();

  bmin[0] = bmi[0];
  bmin[1] = bmi[1];
  bmin[2] = bmi[2];
  bmax[0] = bma[0];
  bmax[1] = bma[1];
  bmax[2] = bma[2];
}



static int is_discrete_ent(struct gmi_model*, struct gmi_ent* e)
{
  (void)e;
  printf("_is_discrete_ent_ not implemented!\n");
  return 0;
}

static void destroy(gmi_model* m)
{
  cap_model* cm = (cap_model*)m;
  free(cm);
}

/* static gmi_model* create_cre(const char* filename) */
/* { */
/*   return gmi_cap_load(filename); */
/* } */

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
  ops.bbox = bbox;
  ops.is_discrete_ent = is_discrete_ent;
  ops.destroy = destroy;
  /* gmi_register(create_cre, "cre"); */
  /* gmi_register(create_native, "xmt_txt"); */
  /* gmi_register(create_native, "x_t"); */
  /* gmi_register(create_native, "sat"); */
}

/* static gmi_model* owned_import(GeometryDatabaseInterface* gi) */
/* { */
/*   gmi_model* m = gmi_import_cap(gi); */
/*   ((cap_model*)m)->owned = true; */
/*   return m; */
/* } */

/* struct gmi_model* gmi_cap_load(const char* creFileName) */
/* { */
/*   if (!gmi_has_ext(creFileName, "cre")) */
/*     gmi_fail("gmi_cap_load: cre file must have .cre extension"); */
/*   const std::string gdbName("Geometry Database : SMLIB");// Switch Create with SMLIB for CAD */
/*   const std::string mdbName("Mesh Database : Create"); */
/*   const std::string adbName("Attribution Database : Create"); */

/*   CapstoneModule  cs("test", gdbName.c_str(), mdbName.c_str(), adbName.c_str()); */

/*   GeometryDatabaseInterface     *g = cs.get_geometry(); */
/*   MeshDatabaseInterface         *m = cs.get_mesh(); */
/*   AppContext                    *c = cs.get_context(); */

/*   PCU_ALWAYS_ASSERT(g); */
/*   PCU_ALWAYS_ASSERT(m); */
/*   PCU_ALWAYS_ASSERT(c); */

/*   v_string filenames; */
/*   filenames.push_back(creFileName); */

/*   M_GModel gmodel = cs.load_files(filenames); */

/*   int numBreps; */
/*   g->get_num_breps(numBreps); */
/*   PCU_ALWAYS_ASSERT(numBreps == 1); */

/*   return owned_import(g); */
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
