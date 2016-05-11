#ifndef APF_BOX_H
#define APF_BOX_H

#include <apfMesh2.h>
#include <gmi_base.h>

namespace apf {

struct Indices
{
  Indices();
  Indices(int a, int b, int c);
  int x, y, z;
  int& operator[](int i);
  Indices operator+(Indices oi);
  Indices operator*(int s);
  static Indices unit(int d);
};

struct Grid {
  Grid(int nx, int ny, int nz);
  int total();
  Indices out(int i);
  int in(Indices is);
  Indices size;
  int stride[4];
};

struct BoxBuilder
{
  Grid grid;
  Grid mgrid;
  int dim;
  double w[3];
  bool is_simplex;
  struct { int dim; int tag; } modelTable[27];
  int modelCounts[4];
  Mesh2* m;
  std::vector<MeshEntity*> v;
  BoxBuilder(int nx, int ny, int nz,
      double wx, double wy, double wz,
      bool is);
  void formModelTable();
  void addModelUse(gmi_base* gb, agm_bdry ab, Indices di);
  gmi_model* buildModel();
  int getModelIndex(int i, int d);
  Indices getModelIndices(Indices vi);
  ModelEntity* getModelEntity(Indices mi);
  void buildCellVert(int i);
  MeshEntity* getVert(Indices vi);
  void buildCellEdges(int i);
  void buildTriangles(MeshEntity* fv[4], ModelEntity* me);
  void buildFace(MeshEntity* fv[4], ModelEntity* me);
  void buildCellFaces(int i);
  void buildTets(MeshEntity* rv[8], ModelEntity* me);
  void buildHex(MeshEntity* rv[8], ModelEntity* me);
  void buildRegion(MeshEntity* rv[8], ModelEntity* me);
  void buildCellRegion(int i);
  void buildCell(int i, int d);
  void buildDimension(int d);
  void buildMeshAndModel();
};

/** \brief create a box from an MDS mesh
  \param nx number of x elements
  \param ny number of y elements
  \param nz number of z elements
  \param wx x dimension width
  \param wy y dimension width
  \param wz z dimension width
  \param is true = simplical mesh, false = quad/hex
  \details set ny,nz=0 for a 1D mesh, set nz=0 for a 2D mesh */
Mesh2* makeMdsBox(
    int nx, int ny, int nz, double wx, double wy, double wz, bool is);

}

#endif
