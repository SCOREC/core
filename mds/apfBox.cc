#include "apfBox.h"

#include "apfMDS.h"
#include <gmi_lookup.h>
#include <pcu_util.h>
#include <lionPrint.h>
#include <cstdlib>

namespace apf {

Indices::Indices() : x(0), y(0), z(0) {}

Indices::Indices(int a, int b, int c): x(a),y(b),z(c) {}

int& Indices::operator[](int i)
{
  switch (i) {
    case 0: return x;
    case 1: return y;
    case 2: return z;
    default: lion_oprint(1,"i must be in {0,1,2}"); abort(); return x;
  }
}

Indices Indices::operator+(Indices oi)
{
  return Indices(x + oi.x, y + oi.y, z + oi.z);
}

Indices Indices::operator*(int s)
{
  return Indices(x * s, y * s, z * s);
}

Indices Indices::unit(int d)
{
  Indices i(0,0,0);
  i[d] = 1;
  return i;
}

Grid::Grid(int nx, int ny, int nz):
  size(nx, ny, nz)
{
  stride[0] = 1;
  for (int i = 0; i < 3; ++i)
    stride[i + 1] = stride[i] * size[i];
}

int Grid::total() {return stride[3];}

Indices Grid::out(int i)
{
  Indices is;
  for (int j = 0; j < 3; ++j)
    is[j] = (i % stride[j + 1]) / stride[j];
  return is;
}

int Grid::in(Indices is)
{
  int i = 0;
  for (int j = 0; j < 3; ++j)
    i += is[j] * stride[j];
  return i;
}

BoxBuilder::BoxBuilder(int nx, int ny, int nz,
      double wx, double wy, double wz,
      bool is, pcu::PCU *PCUObj):
  grid(nx + 1, ny + 1, nz + 1),
  mgrid(nx ? 3 : 1, ny ? 3 : 1, nz ? 3 : 1)
{
  for (dim = 0; dim < 3 && grid.size[dim] > 1; ++dim);
  w[0] = nx ? (wx / nx) : 0;
  w[1] = ny ? (wy / ny) : 0;
  w[2] = nz ? (wz / nz) : 0;
  is_simplex = is;
  formModelTable();
  gmi_model* gm = buildModel();
  m = makeEmptyMdsMesh(gm, dim, false, PCUObj);
  v.resize(grid.total());
  buildMeshAndModel();
}

void BoxBuilder::formModelTable()
{
  int nd[4] = {0,0,0,0};
  for (int i = 0; i < mgrid.total(); ++i) {
    Indices mi = mgrid.out(i);
    int mdim = 0;
    for (int j = 0; j < 3; ++j)
      if (mi[j] == 1)
        ++mdim;
    modelTable[i].dim = mdim;
    modelTable[i].tag = nd[mdim]++;
  }
  for (int i = 0; i < 4; ++i)
    modelCounts[i] = nd[i];
}

int BoxBuilder::getModelIndex(int i, int d)
{
  if (i == 0)
    return 0;
  if (i == grid.size[d] - 1)
    return 2;
  return 1;
}

Indices BoxBuilder::getModelIndices(Indices vi)
{
  Indices mi;
  for (int i = 0; i < 3; ++i)
    mi[i] = getModelIndex(vi[i], i);
  return mi;
}

ModelEntity* BoxBuilder::getModelEntity(Indices mi)
{
  int mj = mgrid.in(mi);
  int mdim = modelTable[mj].dim;
  int mtag = modelTable[mj].tag;
  return m->findModelEntity(mdim, mtag);
}

void BoxBuilder::addModelUse(gmi_base* gb, agm_bdry ab, Indices di)
{
  int ddim = modelTable[mgrid.in(di)].dim;
  int dtag = modelTable[mgrid.in(di)].tag;
  agm_ent ad = gmi_look_up(gb->lookup, agm_type_from_dim(ddim), dtag);
  agm_add_use(gb->topo, ab, ad);
}

gmi_model* BoxBuilder::buildModel()
{
  /* plain malloc because gmi_destroy calls plain free */
  gmi_base* gb = (gmi_base*) malloc(sizeof(*gb));
  gb->model.ops = &gmi_base_ops;
  gmi_base_init(gb);
  agm_ent_type at;
  agm_ent ae;
  agm_bdry ab;
  Indices di;
  for (int i = 0; i < mgrid.total(); ++i) {
    int mdim = modelTable[i].dim;
    int mtag = modelTable[i].tag;
    at = agm_type_from_dim(mdim);
    ae = agm_add_ent(gb->topo, at);
    gmi_set_lookup(gb->lookup, ae, mtag);
  }
  for (int i = 0; i < 4; ++i) {
    at = agm_type_from_dim(i);
    gmi_freeze_lookup(gb->lookup, at);
    gb->model.n[i] = agm_ent_count(gb->topo, at);
    PCU_ALWAYS_ASSERT(gb->model.n[i] == modelCounts[i]);
  }
  for (int i = 0; i < mgrid.total(); ++i) {
    int mdim = modelTable[i].dim;
    int mtag = modelTable[i].tag;
    if (mdim == 0)
      continue;
    at = agm_type_from_dim(mdim);
    ae = gmi_look_up(gb->lookup, at, mtag);
    Indices mi = mgrid.out(i);
    ab = agm_add_bdry(gb->topo, ae);
    for (int j = 0; j < 3; ++j)
      if (mi[j] == 1) {
        di = mi;
        di[j] = 0;
        addModelUse(gb, ab, di);
        di[j] = 2;
        addModelUse(gb, ab, di);
      }
  }
  return &gb->model;
}

void BoxBuilder::buildCellVert(int i)
{
  Indices vi = grid.out(i);
  v[i] = m->createVert(getModelEntity(getModelIndices(vi)));
  Vector3 pt(w[0] * vi.x, w[1] * vi.y, w[2] * vi.z);
  m->setPoint(v[i], 0, pt);
}

MeshEntity* BoxBuilder::getVert(Indices vi)
{
  return v.at(grid.in(vi));
}

void BoxBuilder::buildCellEdges(int i)
{
  Indices vi = grid.out(i);
  Indices mi = getModelIndices(vi);
  MeshEntity* ev[2];
  ev[0] = getVert(vi);
  for (int j = 0; j < 3; ++j) {
    if (mi[j] == mgrid.size[j] - 1)
      continue;
    ev[1] = getVert(vi + Indices::unit(j));
    Indices emi = mi;
    emi[j] = 1;
    ModelEntity* me = getModelEntity(emi);
    buildElement(m, me, Mesh::EDGE, ev);
  }
}

void BoxBuilder::buildTriangles(MeshEntity* fv[4], ModelEntity* me)
{
  MeshEntity* tv[3];
  tv[0] = fv[0]; tv[1] = fv[1]; tv[2] = fv[2];
  buildElement(m, me, Mesh::TRIANGLE, tv);
  tv[0] = fv[2]; tv[1] = fv[3]; tv[2] = fv[0];
  buildElement(m, me, Mesh::TRIANGLE, tv);
}

void BoxBuilder::buildFace(MeshEntity* fv[4], ModelEntity* me)
{
  if (is_simplex)
    buildTriangles(fv, me);
  else
    buildElement(m, me, Mesh::QUAD, fv);
}

void BoxBuilder::buildCellFaces(int i)
{
  Indices vi = grid.out(i);
  Indices mi = getModelIndices(vi);
  MeshEntity* fv[4];
  fv[0] = getVert(vi);
  for (int jx = 0; jx < 3; ++jx) {
    int jy = (jx + 1) % 3;
    if (mi[jx] == mgrid.size[jx] - 1 ||
        mi[jy] == mgrid.size[jy] - 1)
      continue;
    fv[1] = getVert(vi + Indices::unit(jx));
    fv[2] = getVert(vi + Indices::unit(jx) + Indices::unit(jy));
    fv[3] = getVert(vi + Indices::unit(jy));
    Indices emi = mi;
    emi[jx] = 1;
    emi[jy] = 1;
    ModelEntity* me = getModelEntity(emi);
    buildFace(fv, me);
  }
}

void BoxBuilder::buildTets(MeshEntity* rv[8], ModelEntity* me)
{
  static int const tet_verts[6][4] = {
  {0,1,2,6},
  {0,2,3,6},
  {0,3,7,6},
  {0,7,4,6},
  {0,4,5,6},
  {0,5,1,6}};
  for (int i = 0; i < 6; ++i) {
    MeshEntity* tv[4];
    for (int j = 0; j < 4; ++j)
      tv[j] = rv[tet_verts[i][j]];
    buildElement(m, me, Mesh::TET, tv);
  }
}

void BoxBuilder::buildHex(MeshEntity* rv[8], ModelEntity* me)
{
  buildElement(m, me, Mesh::HEX, rv);
}

void BoxBuilder::buildRegion(MeshEntity* rv[8], ModelEntity* me)
{
  if (is_simplex)
    buildTets(rv, me);
  else
    buildHex(rv, me);
}

void BoxBuilder::buildCellRegion(int i)
{
  Indices vi = grid.out(i);
  Indices mi = getModelIndices(vi);
  if (mi.x == mgrid.size.x - 1 ||
      mi.y == mgrid.size.y - 1 ||
      mi.z == mgrid.size.z - 1)
    return;
  MeshEntity* rv[8];
  rv[0] = getVert(Indices(vi.x + 0, vi.y + 0, vi.z + 0));
  rv[1] = getVert(Indices(vi.x + 1, vi.y + 0, vi.z + 0));
  rv[2] = getVert(Indices(vi.x + 1, vi.y + 1, vi.z + 0));
  rv[3] = getVert(Indices(vi.x + 0, vi.y + 1, vi.z + 0));
  rv[4] = getVert(Indices(vi.x + 0, vi.y + 0, vi.z + 1));
  rv[5] = getVert(Indices(vi.x + 1, vi.y + 0, vi.z + 1));
  rv[6] = getVert(Indices(vi.x + 1, vi.y + 1, vi.z + 1));
  rv[7] = getVert(Indices(vi.x + 0, vi.y + 1, vi.z + 1));
  ModelEntity* me = getModelEntity(Indices(1,1,1));
  buildRegion(rv, me);
}

void BoxBuilder::buildCell(int i, int d)
{
  switch (d) {
    case 0:
      buildCellVert(i);
      return;
    case 1:
      buildCellEdges(i);
      return;
    case 2:
      buildCellFaces(i);
      return;
    case 3:
      buildCellRegion(i);
      return;
  };
}

void BoxBuilder::buildDimension(int d)
{
  for (int i = 0; i < grid.total(); ++i)
    buildCell(i, d);
}

void BoxBuilder::buildMeshAndModel()
{
  for (int d = 0; d <= dim; ++d)
    buildDimension(d);
  m->acceptChanges();
}

Mesh2* makeMdsBox(
    int nex, int ney, int nez,
    double wx, double wy, double wz, bool is, pcu::PCU *PCUObj)
{
  BoxBuilder bb(nex, ney, nez, wx, wy, wz, is, PCUObj);
  return bb.m;
}

gmi_model* makeMdsBoxModel(
    int nex, int ney, int nez,
    double wx, double wy, double wz, bool is, pcu::PCU *PCUObj)
{
  BoxBuilder bb(nex, ney, nez, wx, wy, wz, is, PCUObj);
  return bb.buildModel();
}


}
