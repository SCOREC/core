/*
 * Copyright 2015 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "crv.h"
#include "PCU.h"
#include <sstream>
#include <fstream>

namespace crv {

/*
 * Checks if an entity touches a boundary,
 * not just is actually on a boundary
 */
static bool isBoundaryEntity(apf::Mesh* m, apf::MeshEntity* e){
  apf::ModelEntity* g = m->toModel(e);
  int eDim = apf::getDimension(m,e);
  int mDim = m->getDimension();
  int gDim = m->getModelType(g);
  if(gDim != mDim) return true;
  apf::Downward down;
  for(int d = 1; d < eDim; ++d){
    int dDim = m->getDownward(e,d,down);
    for(int i = 0; i < dDim; ++i)
      if(m->getModelType(m->toModel(down[i])) != mDim) return true;
  }
  return false;
}

static int countBoundaryEntities(apf::Mesh* m, int type)
{
  int n = 0;
  apf::MeshIterator* it = m->begin(apf::Mesh::typeDimension[type]);
  apf::MeshEntity* e;
  while ((e = m->iterate(it))) {
    if(m->getType(e) == type && isBoundaryEntity(m,e))
      n++;
  }
  m->end(it);
  return n;
}

static void writePointConnectivity(std::ostream& file, int nPoints)
{
  file << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
  for(int i = 0; i < nPoints; ++i)
    file << i << '\n';
  file << "</DataArray>\n";
}

/*
 * Simple subdivision on an edge into n smaller edges
 */
static void writeEdgeConnectivity(std::ostream& file, apf::Mesh* m, int n)
{
  file << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
  apf::MeshEntity* e;
  apf::MeshIterator* it = m->begin(1);
  int num = 0;
  while ((e = m->iterate(it)))
  {
    if (m->getModelType(m->toModel(e)) != m->getDimension()) {
      for (int i=0; i < n; ++i)
        file << num+i << ' ' << num+i+1 << '\n';
      num += n+1;
    } else {
      file << num << ' ' << num+1 << '\n';
      num += 2;
    }
  }
  m->end(it);
  file << "</DataArray>\n";
}

/*
This will subdivide curved faces by splitting
a triangle into n*n faces

The points are connected from the
last node, onward, for example, for nSplit = 3,
9
8 7
6 5 4
3 2 1 0
where the connectivity in writeFaceConnectivity
will first do the first 6 (n*(n+1)/2)
9 8 7
8 6 5
7 5 4
6 3 2
5 2 1
4 1 0
followed by the last 3 (n*(n-1)/2)
5 8 7
2 6 5
1 5 4
*/
static void writeTriangleConnectivity(std::ostream& file, apf::Mesh* m, int n)
{
  file << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
  apf::MeshEntity* e;
  apf::MeshIterator* it = m->begin(2);
  int num = 0;
  while ((e = m->iterate(it)))
  {
    if(isBoundaryEntity(m,e)) {
      int index = (n+1)*(n+2)/2-1;
      for (int i=0; i < n; ++i){
        for (int j=0; j < i+1; ++j){
          file << num+index-j << ' ' << num+index-j-(i+1)
               << ' ' << num+index-j-(i+2)<< '\n';
        }
        index-=1+i;
      }
      index = (n+1)*(n+2)/2-5;
      for (int i=1; i < n; ++i){
        for (int j=0; j < i; ++j){
          file << num+index-j << ' ' << num+index-j+i+2
               << ' ' << num+index-j+i+1 << '\n';
        }
        index-=2+i;
      }
      num += (n+1)*(n+2)/2;
    } else {
      file << num << ' ' << num+1 << ' ' << num+2 << '\n';
      num += 3;
    }
  }
  m->end(it);
  file << "</DataArray>\n";
}

/*
 * Tets are subdivided into triangular prisms and tets
 * by choosing one face, subdividing it into triangles as above
 * and then creating layers until the opposite vertex is approached
 * which is connected to that vertex with tets
 */
static void writeTetConnectivity(std::ostream& file, apf::Mesh* m, int n)
{
  file << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
  apf::MeshEntity* e;
  apf::MeshIterator* it = m->begin(3);
  int num = 0, apex = n*(n+1)*(n+2)/2;
  while ((e = m->iterate(it)))
  {
    if(isBoundaryEntity(m,e)) {
      for(int layer = 0; layer < (n-1); ++layer){
        int index = (n+1)*(n+2)/2-1;
        for (int i=0; i < n; ++i){
          for (int j=0; j < i+1; ++j){
            file << (num+index-j) << ' '
                 << (num+index-j-(i+1)) << ' '
                 << (num+index-j-(i+2)) << ' '
                 << (n+1)*(n+2)/2+(num+index-j) << ' '
                 << (n+1)*(n+2)/2+(num+index-j-(i+1)) << ' '
                 << (n+1)*(n+2)/2+(num+index-j-(i+2)) << '\n';
          }
          index-=1+i;
        }
        index = (n+1)*(n+2)/2-5;
        for (int i=1; i < n; ++i){
          for (int j=0; j < i; ++j){
            file << (num+index-j) << ' '
                 << (num+index-j+i+2) << ' '
                 << (num+index-j+i+1) << ' '
                 << (n+1)*(n+2)/2+(num+index-j) << ' '
                 << (n+1)*(n+2)/2+(num+index-j+i+2) << ' '
                 << (n+1)*(n+2)/2+(num+index-j+i+1) << '\n';
          }
          index-=2+i;
        }
        num += (n+1)*(n+2)/2;
      }
      num += 1;
      int index = (n+1)*(n+2)/2-2;
      for (int i=0; i < n; ++i){
        for (int j=0; j < i+1; ++j){
          file << (num+index-j) << ' '
              << (num+index-j-(i+1)) << ' '
              << (num+index-j-(i+2)) << ' '
              << apex << '\n';
        }
        index-=1+i;
      }
      index = (n+1)*(n+2)/2-6;
      for (int i=1; i < n; ++i){
        for (int j=0; j < i; ++j){
          file << (num+index-j) << ' '
              << (num+index-j+i+2) << ' '
              << (num+index-j+i+1) << ' '
              << apex << '\n';
        }
        index-=2+i;
      }
      num += (n+1)*(n+2)/2;
      apex += n*(n+1)*(n+2)/2+1;
    } else {
      file << num << ' ' << num+1 << ' ' << num+2 << ' ' << num+3 << '\n';
      num += 4;
    }
  }
  m->end(it);
  file << "</DataArray>\n";
}

static void writeConnectivity(std::ostream& file, apf::Mesh* m, int type, int n)
{
  switch (type) {
    case apf::Mesh::VERTEX:
      writePointConnectivity(file,n);
      break;
    case apf::Mesh::EDGE:
      writeEdgeConnectivity(file,m,n);
      break;
    case apf::Mesh::TRIANGLE:
      writeTriangleConnectivity(file,m,n);
      break;
    case apf::Mesh::TET:
      writeTetConnectivity(file,m,n);
      break;
    default:
      fail("can only write curved VTU files for control points, "
          "edges, triangles, and tets");
      break;
  }
}

static void writeOffsets(std::ostream& file, int type, int nCells)
{
  file << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
  int o = 0;
  for (int i=0; i < nCells; ++i){
    o += apf::Mesh::typeDimension[type]+1;
    file << o << '\n';
  }
  file << "</DataArray>\n";
}

static void writeTypes(std::ostream& file, int type, int nCells)
{
  file << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
  static int vtkTypes[3] = {1,3,5};
  for (int i=0; i < nCells; ++i)
    file << vtkTypes[type] << '\n';
  file << "</DataArray>\n";
}

static void writePoint(std::ostream& file, apf::Vector3 & pt)
{
  for (int j=0; j < 3; ++j)
    file << pt[j] << ' ';
  file << '\n';
}

static void writeStart(std::ostream& file, int nPoints, int nCells)
{
  file << "<VTKFile type=\"UnstructuredGrid\">\n";
  file << "<UnstructuredGrid>\n";
  file << "<Piece NumberOfPoints=\"" << nPoints;
  file << "\" NumberOfCells=\"" << nCells;
  file << "\">\n";
}

static void writeCells(std::ostream& file, apf::Mesh* m,
    int type, int nSplit, int nCells)
{
  file << "<Cells>\n";
  writeConnectivity(file,m,type,nSplit);
  writeOffsets(file,type,nCells);
  writeTypes(file,type,nCells);
  file << "</Cells>\n";
}

static void writeEnd(std::ostream& file)
{
  file << "</Piece>\n";
  file << "</UnstructuredGrid>\n";
  file << "</VTKFile>\n";
}

static void writeTetJacobianData(std::ostream& file, apf::Mesh* m, int n)
{
  file << "<PointData>\n";
  file << "<DataArray type=\"Float64\" Name=\"detJacobian\" "
       << "NumberOfComponents=\"1\" format=\"ascii\">\n";

  apf::MeshIterator* it = m->begin(3);
  apf::MeshEntity* e;
  apf::MeshEntity* down[4];
  apf::Vector3 p,pt;
  apf::Matrix3x3 J;
  apf::Vector3 params[4] = {apf::Vector3(0,0,0),
        apf::Vector3(1,0,0),apf::Vector3(0,1,0),apf::Vector3(0,0,1)};

  while ((e = m->iterate(it))) {
    apf::MeshElement* me = apf::createMeshElement(m,e);
    if(isBoundaryEntity(m,e)){
      for (int k = 0; k < n; ++k){
        p[2] = 1.*k/n;
        for (int j = 0; j <= n; ++j){
          p[1] = 1.*j/n*(1.-p[2]);
          for (int i = 0; i <= n-j; ++i){
            p[0] = 1.*i/n*(1.-p[2]);
            apf::getJacobian(me,p,J);
            file << apf::getDeterminant(J) << '\n';
          }
        }
      }
      p = apf::Vector3(0.,0.,1.); // apex point
      apf::getJacobian(me,p,J);
      file << apf::getDeterminant(J) << '\n';
    } else {
      m->getDownward(e,0,down);
      for(int i = 0; i < 4; ++i){
        apf::getJacobian(me,params[i],J);
        file << apf::getDeterminant(J) << '\n';
      }
    }
    apf::destroyMeshElement(me);
  }

  file << "</DataArray>\n";
  file << "</PointData>\n";
}

static void writeEdgeVtuFiles(apf::Mesh* m, int n, const char* prefix)
{
  std::stringstream ss;
  ss << prefix << PCU_Comm_Self() << "_"
     << m->getShape()->getOrder()
      << "_edges.vtu";

  int nBoundaryEnts = countBoundaryEntities(m,apf::Mesh::EDGE);

  int nPoints = 2*m->count(1) + nBoundaryEnts*(n-1);
  int nCells = m->count(1) + nBoundaryEnts*(n-1);

  std::string fileName = ss.str();
  std::ofstream file(fileName.c_str());
  assert(file.is_open());

  writeStart(file,nPoints,nCells);

  file << "<Points>\n";
  file << "<DataArray type=\"Float64\" Name=\"coordinates\" "
      "NumberOfComponents=\"3\" format=\"ascii\">\n";

  apf::MeshIterator* it = m->begin(1);
  apf::MeshEntity* e;
  apf::MeshEntity* v[2];
  apf::Vector3 p,pt;
  while ((e = m->iterate(it))) {
    if (isBoundaryEntity(m,e)){
      apf::Element* elem =
          apf::createElement(m->getCoordinateField(),e);
      for (int i = 0; i <= n; ++i){
        p[0] = 2.*i/n-1.;
        apf::getVector(elem,p,pt);
        writePoint(file,pt);
      }
      apf::destroyElement(elem);

    } else {
      m->getDownward(e,0,v);
      for(int i = 0; i < 2; ++i){
        m->getPoint(v[i],0,pt);
        writePoint(file,pt);
      }
    }
  }
  m->end(it);
  file << "</DataArray>\n";
  file << "</Points>\n";
  writeCells(file,m,apf::Mesh::EDGE,n,nCells);
  writeEnd(file);
}

static void writeTriangleVtuFiles(apf::Mesh* m, int n, const char* prefix)
{
  std::stringstream ss;
  ss << prefix << PCU_Comm_Self() << "_"
     << m->getShape()->getOrder()
     << "_faces.vtu";

  int nBoundaryEnts = countBoundaryEntities(m,apf::Mesh::TRIANGLE);

  int nPoints = 3*m->count(2) + nBoundaryEnts*((n+1)*(n+2)/2-3);
  int nCells = m->count(2) + nBoundaryEnts*(n*n-1);

  std::string fileName = ss.str();
  std::ofstream file(fileName.c_str());
  assert(file.is_open());

  writeStart(file,nPoints,nCells);

  file << "<Points>\n";
  file << "<DataArray type=\"Float64\" Name=\"coordinates\" "
      "NumberOfComponents=\"3\" format=\"ascii\">\n";

  apf::MeshIterator* it = m->begin(2);
  apf::MeshEntity* e;
  apf::MeshEntity* down[3];
  apf::Vector3 p,pt;
  while ((e = m->iterate(it))) {
    if(isBoundaryEntity(m,e)){
      apf::Element* elem =
          apf::createElement(m->getCoordinateField(),e);
      for (int j = 0; j <= n; ++j){
        p[1] = 1.*j/n;
        for (int i = 0; i <= n-j; ++i){
          p[0] = 1.*i/n;
          apf::getVector(elem,p,pt);
          writePoint(file,pt);
        }
      }
      apf::destroyElement(elem);
    } else {
      m->getDownward(e,0,down);
      for(int i = 0; i < 3; ++i){
        m->getPoint(down[i],0,pt);
        writePoint(file,pt);
      }
    }
  }
  m->end(it);
  file << "</DataArray>\n";
  file << "</Points>\n";
  writeCells(file,m,apf::Mesh::TRIANGLE,n,nCells);
  writeEnd(file);
}

static void writeTetVtuFiles(apf::Mesh* m, int n, const char* prefix)
{
  std::stringstream ss;
  ss << prefix << PCU_Comm_Self() << "_"
     << m->getShape()->getOrder()
     << "_tets.vtu";

  int nBoundaryEnts = countBoundaryEntities(m,apf::Mesh::TET);

  int nPoints = 4*m->count(3) + nBoundaryEnts*(n*(n+1)*(n+2)/2+1-4);
  int nCells = m->count(3) + nBoundaryEnts*(n*n*n-1);

  std::string fileName = ss.str();
  std::ofstream file(fileName.c_str());
  assert(file.is_open());

  writeStart(file,nPoints,nCells);
  file << "<Points>\n";
  file << "<DataArray type=\"Float64\" Name=\"coordinates\" "
      "NumberOfComponents=\"3\" format=\"ascii\">\n";

  apf::MeshIterator* it = m->begin(3);
  apf::MeshEntity* e;
  apf::MeshEntity* down[4];
  apf::Vector3 p,pt;
  while ((e = m->iterate(it))) {
    if(isBoundaryEntity(m,e)){
      apf::Element* elem =
          apf::createElement(m->getCoordinateField(),e);
      for (int k = 0; k < n; ++k){
        p[2] = 1.*k/n;
        for (int j = 0; j <= n; ++j){
          p[1] = 1.*j/n*(1.-p[2]);
          for (int i = 0; i <= n-j; ++i){
            p[0] = 1.*i/n*(1.-p[2]);
            apf::getVector(elem,p,pt);
            writePoint(file,pt);
          }
        }
      }
      p = apf::Vector3(0.,0.,1.); // apex point
      apf::getVector(elem,p,pt);
      writePoint(file,pt);
      apf::destroyElement(elem);
    } else {
      m->getDownward(e,0,down);
      for(int i = 0; i < 4; ++i){
        m->getPoint(down[i],0,pt);
        writePoint(file,pt);
      }
    }
  }
  file << "</DataArray>\n";
  file << "</Points>\n";
  file << "<Cells>\n";

  writeTetConnectivity(file,m,n);

  file << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
  int o = 0;
  it = m->begin(3);
  while ((e = m->iterate(it))) {
    if(isBoundaryEntity(m,e)){
      for(int i = 0; i < n*n*(n-1); ++i){
        o += 6;
        file << o << '\n';
      }
      for(int i = 0; i < n*n; ++i){
        o += 4;
        file << o << '\n';
      }
    } else {
      o += 4;
      file << o << '\n';
    }
  }
  file << "</DataArray>\n";

  file << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
  it = m->begin(3);
  while ((e = m->iterate(it))) {
    if(isBoundaryEntity(m,e)){
      for(int i = 0; i < n*n*(n-1); ++i)
        file << 13 << '\n';

      for(int i = 0; i < n*n; ++i)
        file << 10 << '\n';

    } else
      file << 10 << '\n';
  }

  file << "</DataArray>\n";
  file << "</Cells>\n";
  writeTetJacobianData(file,m,n);
  writeEnd(file);
  m->end(it);
}

static void writeControlPointVtuFiles(apf::Mesh* m, const char* prefix)
{
  std::stringstream ss;
  ss << prefix << PCU_Comm_Self() << "_"
     << m->getShape()->getOrder()
     << "_controlPoints.vtu";
  int nPoints = 0;
  for (int t = 0; t < apf::Mesh::TYPES; ++t)
    nPoints += m->getShape()->countNodesOn(t)*countBoundaryEntities(m,t);

  std::string fileName = ss.str();
  std::ofstream file(fileName.c_str());
  assert(file.is_open());

  writeStart(file,nPoints,nPoints);
  file << "<Points>\n";
  file << "<DataArray type=\"Float64\" Name=\"coordinates\" "
      "NumberOfComponents=\"3\" format=\"ascii\">\n";

  for (int t = 0; t < apf::Mesh::TYPES; ++t){
    apf::MeshIterator* it = m->begin(apf::Mesh::typeDimension[t]);
    apf::MeshEntity* e;
    apf::Vector3 pt;
    while ((e = m->iterate(it))) {
      if(isBoundaryEntity(m,e)){
        for(int i = 0; i < m->getShape()->countNodesOn(t); ++i){
          m->getPoint(e,i,pt);
          writePoint(file,pt);
        }
      }
    }
    m->end(it);
  }
  file << "</DataArray>\n";
  file << "</Points>\n";
  writeCells(file,m,apf::Mesh::VERTEX,nPoints,nPoints);
  writeEnd(file);
}

void writeCurvedVtuFiles(apf::Mesh* m, int type, int n, const char* prefix)
{
  switch (type) {
    case apf::Mesh::VERTEX:
      writeControlPointVtuFiles(m,prefix);
      break;
    case apf::Mesh::EDGE:
      writeEdgeVtuFiles(m,n,prefix);
      break;
    case apf::Mesh::TRIANGLE:
      writeTriangleVtuFiles(m,n,prefix);
      break;
    case apf::Mesh::TET:
      writeTetVtuFiles(m,n,prefix);
      break;
    default:
      break;
  }
}

} //namespace crv
