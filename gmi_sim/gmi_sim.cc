/****************************************************************************** 

  Copyright 2014 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/
#include <PCU.h>
#include "gmi_sim.h"
#include <gmi.h>
#include <stdlib.h>
#include <gmi.h>
#include <SimModel.h>
#include <vector>

#include "gmi_sim_config.h"

#ifdef SIM_PARASOLID
#include "SimParasolidKrnl.h"
#endif

#ifdef SIM_ACIS
#include "SimAcisKrnl.h"
#endif

struct sim_model {
  struct gmi_model model;
  SGModel* sim;
  bool owned;
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


/* apparently much of what Simmetrix does
   uses a custom memory allocator which is not
   thread-safe. This is a relic of the old scorec
   codebase that simmetrix built on.
   hopefully they will get rid of it..
   until then we protect these code segments
   with a spinlock */

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

static gmi_ent* next(gmi_model*, gmi_iter* i)
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

static void end(gmi_model*, gmi_iter* i)
{
  sim_iter* si;
  si = (sim_iter*)i;
  if (si->dim == 0)
    GVIter_delete(si->i.v);
  else if (si->dim == 1)
    GEIter_delete(si->i.e);
  else if (si->dim == 2)
    GFIter_delete(si->i.f);
  else if (si->dim == 3)
    GRIter_delete(si->i.r);
  free(si);
}

static int get_dim(gmi_model*, gmi_ent* e)
{
  return GEN_type((pGEntity)e);
}

static int get_tag(gmi_model*, gmi_ent* e)
{
  return GEN_tag((pGEntity)e);
}

static gmi_ent* find(gmi_model* m, int dim, int tag)
{
  sim_model* mm = (sim_model*)m;
  return (gmi_ent*)GM_entityByTag(mm->sim, dim, tag);
}

static gmi_set* plist_to_set(pPList l)
{
  gmi_set* s = gmi_make_set(PList_size(l));
  for (int i = 0; i < s->n; ++i)
    s->e[i] = (gmi_ent*)PList_item(l, i);
  PList_delete(l);
  return s;
}

/* GF_regions removes duplicates (in some versions of the code)
   which destroys valuable non-manifold information for us.
   So, custom code here */
static int face_regions2(pGFace face, pGRegion regions[2])
{
  int nregions = 0;
  for (int i = 0; i < 2; ++i) {
    pGFaceUse use = GF_use(face, i);
    if (use) {
      pGRegion region = GFU_region(use);
      if (region)
        regions[nregions++] = region;
    }
  }
  return nregions;
}

static gmi_set* face_regions(pGFace face)
{
  pGRegion regions[2];
  int nregions = face_regions2(face, regions);
  gmi_set* s = gmi_make_set(nregions);
  for (int i = 0; i < nregions; ++i)
    s->e[i] = (gmi_ent*)regions[i];
  return s;
}

/* same story as above, want to preserve duplicates
   for non-manifold faces.
 unfortunately, simmetrix doesn't expose a shell api,
 so we have to do this horrible "check whether each
 face should actually have shown up twice" thing. */
static gmi_set* region_faces(pGRegion region)
{
  pPList unique_faces;
  std::vector<pGFace> faces;
  unique_faces = GR_faces(region);
  faces.reserve(PList_size(unique_faces));
  for (int i = 0; i < PList_size(unique_faces); ++i) {
    pGFace face = (pGFace)PList_item(unique_faces, i);
    faces.push_back(face);
    pGRegion regions[2];
    int nregions = face_regions2(face, regions);
    if ((nregions == 2) && (regions[0] == regions[1]))
      faces.push_back(face);
  }
  PList_delete(unique_faces);
  gmi_set* s = gmi_make_set(faces.size());
  for (int i = 0; i < s->n; ++i)
    s->e[i] = (gmi_ent*)faces[i];
  return s;
}

/* getting the region adj to an edge. This version
 * does not support non-manifold models
 */
static gmi_set* edge_regions(pGEdge e)
{
  pPList list = GE_faces((pGEdge)e);

  // do the first one outside of the loop
  gmi_set* regions_set = face_regions((pGFace)PList_item(list, 0));
  if (regions_set->n != 1)
    gmi_fail("no support for non-manifold surfaces!\n");
  pGRegion r = (pGRegion)regions_set->e[0];

  // do the rest inside of the loop
  for (int i = 1; i < PList_size(list); i++) {
    regions_set = face_regions((pGFace)PList_item(list, i));
    if (regions_set->n != 1)
      gmi_fail("no support for non-manifold surfaces!\n");
    if (r != (pGRegion)regions_set->e[0])
      gmi_fail("no support for non-manifold surfaces!\n");
  }
  return regions_set;
}

static gmi_set* adjacent(gmi_model* m, gmi_ent* e, int dim)
{
  int edim = gmi_dim(m, e);
  if (edim == 0 && dim == 1)
    return plist_to_set(GV_edges((pGVertex)e));
  if (edim == 1 && dim == 0)
    return plist_to_set(GE_vertices((pGEdge)e));
  if (edim == 1 && dim == 2)
    return plist_to_set(GE_faces((pGEdge)e));
  if (edim == 2 && dim == 0)
    return plist_to_set(GF_vertices((pGFace)e));
  if (edim == 2 && dim == 1)
    return plist_to_set(GF_edges((pGFace)e));
  if (edim == 2 && dim == 3)
    return face_regions((pGFace)e);
  if (edim == 1 && dim == 3)
    return edge_regions((pGEdge)e);
  if (edim == 3 && dim == 2)
    return region_faces((pGRegion)e);
  if (edim == 3 && dim == 4) /* sometimes people just keep looking up */
    return gmi_make_set(0);
  gmi_fail("requested adjacency not available\n");
  return 0;
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
    return;
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

static void closest_point(struct gmi_model* m, struct gmi_ent* e,
    double const from[3], double to[3], double to_p[2])
{
  int md = gmi_dim(m, e);
  if (md == 2)
    GF_closestPoint((pGFace)e,&from[0],&to[0],&to_p[0]);
  else if (md == 1)
    GE_closestPoint((pGEdge)e,&from[0],&to[0],&to_p[0]);
}

static void normal(struct gmi_model* m, struct gmi_ent* e,
    double const p[2], double n[3])
{
  int md = gmi_dim(m, e);
  if (md == 2)
    GF_normal((pGFace)e,&p[0],&n[0]);

}

static void first_derivative(struct gmi_model* m, struct gmi_ent* e,
    double const p[2], double t0[3], double t1[3])
{
  int md = gmi_dim(m, e);
  if (md == 2)
    GF_firstDerivative((pGFace)e,&p[0],&t0[0],&t1[0]);
  else if (md == 1)
    GE_firstDerivative((pGEdge)e,p[0],&t0[0]);
}

static int is_point_in_region(struct gmi_model* m, struct gmi_ent* e,
    double point[3])
{
  gmi_dim(m, e);
  int res = GR_containsPoint((pGRegion)e, &point[0]);
  return res;
}

static int is_in_closure_of(struct gmi_model* m, struct gmi_ent* e,
    struct gmi_ent* et)
{
  int etd = gmi_dim(m, et);
  if (etd == 3)
    return GR_inClosure((pGRegion)et, (pGEntity)e);
  if (etd == 2)
    return GF_inClosure((pGFace)et, (pGEntity)e);
  if (etd == 1)
    return GE_inClosure((pGEdge)et, (pGEntity)e);
  gmi_fail("requested operation is not possible\n");
  return 0;
}

static void destroy(gmi_model* m)
{
  sim_model* mm = (sim_model*)m;
  if (mm->owned)
    GM_release(mm->sim);
  free(mm);
}

static struct gmi_model_ops ops;

static gmi_model* create_smd(const char* filename)
{
  return gmi_sim_load(NULL, filename);
}

static gmi_model* create_native(const char* filename)
{
  return gmi_sim_load(filename, NULL);
}

} //extern "C"

void gmi_sim_start(void)
{
  SimModel_start();
#ifdef SIM_PARASOLID
  SimParasolid_start(1);
#endif
#ifdef SIM_ACIS
  SimAcis_start(1);
#endif
}

void gmi_sim_stop(void)
{
#ifdef SIM_ACIS
  SimAcis_stop(1);
#endif
#ifdef SIM_PARASOLID
  SimParasolid_stop(1);
#endif
  SimModel_stop();
}

void gmi_register_sim(void)
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
  ops.destroy = destroy;
  gmi_register(create_smd, "smd");
  gmi_register(create_native, "xmt_txt");
  gmi_register(create_native, "x_t");
  gmi_register(create_native, "sat");
}

static gmi_model* owned_import(pGModel sm)
{
  gmi_model* m = gmi_import_sim(sm);
  ((sim_model*)m)->owned = true;
  return m;
}

#ifdef SIM_PARASOLID
static pNativeModel load_parasolid(const char* filename)
{
  enum { TEXT_FORMAT = 0 };
  return ParasolidNM_createFromFile(filename, TEXT_FORMAT);
}
#else
static pNativeModel load_parasolid(const char* filename)
{
  (void)filename;
  gmi_fail("recompile with -DSIM_PARASOLID=ON");
}
#endif

#ifdef SIM_ACIS
static pNativeModel load_acis(const char* filename)
{
  enum { TEXT_FORMAT = 0 };
  return AcisNM_createFromFile(filename, TEXT_FORMAT);
}
#else
static pNativeModel load_acis(const char* /*filename*/)
{
  gmi_fail("recompile with -DSIM_ACIS=ON");
}
#endif

struct gmi_model* gmi_sim_load(const char* nativefile, const char* smdfile)
{
  pNativeModel nm;
  if (!nativefile)
    nm = 0;
  else if (gmi_has_ext(nativefile, "sat"))
    nm = load_acis(nativefile);
  else if (gmi_has_ext(nativefile, "xmt_txt"))
    nm = load_parasolid(nativefile);
  else if (gmi_has_ext(nativefile, "x_t"))
    nm = load_parasolid(nativefile);
  else
    gmi_fail("gmi_sim_load: nativefile has bad extension");
  pGModel sm;
  if (!smdfile) {
    if (NM_isAssemblyModel(nm)) {
      pGAModel am = GAM_createFromNativeModel(nm, NULL);
      NM_release(nm);
      sm = GM_createFromAssemblyModel(am, NULL, NULL);
      GM_release(am);
      nm = GM_nativeModel(sm);
    } else
      sm = GM_createFromNativeModel(nm, NULL);
  } else if (gmi_has_ext(smdfile, "smd"))
    sm = GM_load(smdfile, nm, NULL);
  else
    gmi_fail("gmi_sim_load: smdfile has bad extension");
  if (nm)
    NM_release(nm);
  return owned_import(sm);
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
  m->owned = false;
  return &m->model;
}

SGModel* gmi_export_sim(gmi_model* m)
{
  sim_model* mm = (sim_model*)m;
  return mm->sim;
}
