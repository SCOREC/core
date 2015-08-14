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

static int countOwnedEntitiesOfType(apf::Mesh* m, int type)
{
  apf::MeshIterator* it = m->begin(apf::Mesh::typeDimension[type]);
  apf::MeshEntity* e;
  int count = 0;
  while ((e = m->iterate(it)))
    if (m->getType(e)==type && m->isOwned(e))
      ++count;
  m->end(it);
  return count;
}

static std::string getSuffix(int type)
{
  std::stringstream ss;
  switch (type) {
  // control points
  case apf::Mesh::VERTEX:
    ss << "_ctrlPts";
    break;
  case apf::Mesh::EDGE:
    ss << "_edge";
    break;
  case apf::Mesh::TRIANGLE:
    ss << "_tri";
    break;
  case apf::Mesh::TET:
    ss << "_tet";
    break;
  default:
    break;
  }
  return ss.str();
}

static void writePointConnectivity(std::ostream& file, int n)
{
  for(int i = 0; i < n; ++i)
    file << i << '\n';
}

/*
 * Simple subdivision on an edge into n smaller edges
 */
static void writeEdgeConnectivity(std::ostream& file, int c, int n)
{
  int num = 0;
  for(int j = 0; j < c; ++j){
    for(int i = 0; i < n; ++i)
      file << num+i << ' ' << num+i+1 << '\n';
    num += n+1;
  }
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
static void writeTriangleConnectivity(std::ostream& file, int c, int n)
{
  int num = 0;
  for(int k = 0; k < c; ++k){
    int index = (n+1)*(n+2)/2-1;
    for(int i = 0; i < n; ++i){
      for(int j = 0; j < i+1; ++j){
        file << num+index-j << ' ' << num+index-j-(i+1)
             << ' ' << num+index-j-(i+2)<< '\n';
      }
      index -= 1+i;
    }
    index = (n+1)*(n+2)/2-5;
    for (int i = 1; i < n; ++i){
      for (int j = 0; j < i; ++j){
        file << num+index-j << ' ' << num+index-j+i+2
             << ' ' << num+index-j+i+1 << '\n';
      }
      index -= 2+i;
    }
    num += (n+1)*(n+2)/2;
  }
}

/*
 * Tets are subdivided into four hexes, which are then split into more
 * hexes. This gives a more uniform subdivision
 */
static void writeTetConnectivity(std::ostream& file, int c, int n)
{
  int num = 0;
  for(int l = 0; l < c; ++l){
    for(int h = 0; h < 4; ++h){
      for(int k = 0; k < n; ++k){
        for(int j = 0; j < n; ++j){
          for(int i = 0; i < n; ++i){
            file << num+i<< ' ' << num+i+1 << ' ' << num+i+1+n+1 << ' '
                 << num+i+n+1 << ' ' << num+i+(n+1)*(n+1)<< ' '
                 << num+i+1+(n+1)*(n+1) << ' ' << num+i+1+n+1+(n+1)*(n+1)
                 << ' ' << num+i+n+1+(n+1)*(n+1) << '\n';
          }
          num+=n+1;
        }
        num+=n+1;
      }
      num+=(n+1)*(n+1);
    }
  }
}

static void writeConnectivity(std::ostream& file, int type, int c, int n)
{
  file << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
  switch (type) {
    case apf::Mesh::VERTEX:
      writePointConnectivity(file,n);
      break;
    case apf::Mesh::EDGE:
      writeEdgeConnectivity(file,c,n);
      break;
    case apf::Mesh::TRIANGLE:
      writeTriangleConnectivity(file,c,n);
      break;
    case apf::Mesh::TET:
      writeTetConnectivity(file,c,n);
      break;
    default:
      fail("can only write curved VTU files for control points, \
           edges, triangles, and tets");
      break;
  }
  file << "</DataArray>\n";
}

static void writeOffsets(std::ostream& file, int type, int nCells)
{
  file << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
  static int vtkOffsets[apf::Mesh::TYPES] = {
    1,  //parent vertex
    2,  //parent edge
    3,  //parent triangle
    -1, //parent quad
    8, //parent tet, split into hexes, use hex type
    -1,
    -1,
    -1
  };
  int o = 0;
  for (int i=0; i < nCells; ++i){
    o += vtkOffsets[type];
    file << o << '\n';
  }
  file << "</DataArray>\n";
}

static void writeTypes(std::ostream& file, int type, int nCells)
{
  file << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
  assert(type >= 0);
  assert(type < apf::Mesh::TYPES);
  static int vtkTypes[apf::Mesh::TYPES] = {
    1,  //parent vertex
    3,  //parent edge
    5,  //parent triangle
    -1, //parent quad
    12, //parent tet, split into hexes, use hex type
    -1,
    -1,
    -1
  };
  assert(vtkTypes[type] != -1);
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

static void writeCells(std::ostream& file,
    int type, int nEntities, int nSplit, int nCells)
{
  file << "<Cells>\n";
  writeConnectivity(file,type,nEntities,nSplit);
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

static void writeTriJacobianData(std::ostream& file, apf::Mesh* m, int n)
{
  file << "<PointData>\n";
  file << "<DataArray type=\"Float64\" Name=\"detJacobian\" "
       << "NumberOfComponents=\"1\" format=\"ascii\">\n";

  apf::MeshIterator* it = m->begin(2);
  apf::MeshEntity* e;
  apf::Vector3 p;

  apf::Matrix3x3 J;

  while ((e = m->iterate(it))) {
    if(!m->isOwned(e)) continue;
    apf::MeshElement* me =
        apf::createMeshElement(m,e);
    for (int j = 0; j <= n; ++j){
      p[1] = 1.*j/n;
      for (int i = 0; i <= n-j; ++i){
        p[0] = 1.*i/n;
        apf::getJacobian(me,p,J);
        double detJ = apf::getJacobianDeterminant(J,2);
        file << detJ << '\n';
        if(detJ < 0.){
          apf::Vector3 pt;
          apf::getVector((apf::Element*)me,p,pt);
          std::stringstream ss;
          ss << "warning: Tri Jacobian Determinant is negative,  " << detJ
              << " at " << pt << '\n';
          std::string s = ss.str();
          fprintf(stderr, "%s", s.c_str());
        }
      }
    }
    apf::destroyMeshElement(me);
  }
  m->end(it);

  file << "</DataArray>\n";
  file << "</PointData>\n";
}

static void writeTetJacobianData(std::ostream& file, apf::Mesh* m, int n)
{
  file << "<PointData>\n";
  file << "<DataArray type=\"Float64\" Name=\"detJacobian\" "
       << "NumberOfComponents=\"1\" format=\"ascii\">\n";

  apf::MeshIterator* it = m->begin(3);
  apf::MeshEntity* e;
  apf::Vector3 xi,p,pt;
  apf::Matrix3x3 J;
  // first initializing with end points
  apf::Vector3 params[15] = {apf::Vector3(0,0,0),apf::Vector3(1,0,0),
      apf::Vector3(0,1,0),apf::Vector3(0,0,1),apf::Vector3(0,0,0),
      apf::Vector3(0,0,0),apf::Vector3(0,0,0),apf::Vector3(0,0,0),
      apf::Vector3(0,0,0),apf::Vector3(0,0,0),apf::Vector3(0,0,0),
      apf::Vector3(0,0,0),apf::Vector3(0,0,0),apf::Vector3(0,0,0),
      apf::Vector3(0.25,0.25,0.25)};

  for(int i = 0; i < 6; ++i)
    params[i+4] = params[apf::tet_edge_verts[i][0]]*0.5
      + params[apf::tet_edge_verts[i][1]]*0.5;
  for(int i = 0; i < 4; ++i)
    params[i+10] = params[apf::tet_tri_verts[i][0]]*1./3.
      + params[apf::tet_tri_verts[i][1]]*1./3.
      + params[apf::tet_tri_verts[i][2]]*1./3.;

  int hex[4][8] = {{0,4,10,6,7,11,14,13},{1,5,10,4,8,12,14,11},
      {2,6,10,5,9,13,14,12},{9,13,14,12,3,7,11,8}};

  apf::NewArray<double> values;
  apf::EntityShape* shape = apf::getLagrange(1)->getEntityShape(apf::Mesh::HEX);

  while ((e = m->iterate(it))) {
    if(!m->isOwned(e)) continue;
    apf::MeshElement* me = apf::createMeshElement(m,e);
    for(int h = 0; h < 4; ++h){
      for (int k = 0; k <= n; ++k){
        xi[2] = 2.*k/n - 1.;
        for (int j = 0; j <= n; ++j){
          xi[1] = 2.*j/n - 1.;
          for (int i = 0; i <= n; ++i){
            xi[0] = 2.*i/n - 1.;
            shape->getValues(0, 0, xi, values);
            p.zero();
            for(int l = 0; l < 8; ++l)
              p += params[hex[h][l]]*values[l];

            apf::getJacobian(me,p,J);
            double detJ = apf::getDeterminant(J);
            if(detJ < 0.){
              apf::Vector3 pt;
              apf::getVector((apf::Element*)me,p,pt);
              std::stringstream ss;
              ss << "warning: Tet Jacobian Determinant is negative,  " << detJ
                 << " at " << pt << "for " << m->getShape()->getOrder()
                 << " order method\n";
              std::string s = ss.str();
              fprintf(stderr, "%s", s.c_str());
            }
            file << detJ << '\n';
          }
        }
      }
    }
    apf::destroyMeshElement(me);
  }
  m->end(it);

  file << "</DataArray>\n";
  file << "</PointData>\n";
}

static void writePvtuFile(const char* prefix, apf::Mesh* m, int type)
{
  std::stringstream ss;
  ss << prefix << "_" << m->getShape()->getOrder()
     << getSuffix(type) << ".pvtu";
  std::string fileName = ss.str();
  std::ofstream file(fileName.c_str());
  assert(file.is_open());
  file << "<VTKFile type=\"PUnstructuredGrid\">\n";
  file << "<PUnstructuredGrid GhostLevel=\"0\">\n";
  file << "<PPoints>\n";
  file << "<PDataArray type=\"Float64\" Name=\"coordinates\" "
       << "NumberOfComponents=\"3\" format=\"ascii\"/>\n";
  file << "</PPoints>\n";
  file << "<PPointData>\n";
  file << "<PDataArray type=\"Float64\" Name=\"detJacobian\" "
       << "NumberOfComponents=\"1\" format=\"ascii\"/>\n";
  file << "</PPointData>\n";
  for (int i=0; i < PCU_Comm_Peers(); ++i)
  {
    std::stringstream ssPCU;
    ssPCU << prefix << i << "_"
       << m->getShape()->getOrder()
       << getSuffix(type) << ".vtu";
    file << "<Piece Source=\"" << ssPCU.str() << "\"/>\n";
  }

  file << "</PUnstructuredGrid>\n";
  file << "</VTKFile>\n";
}

static void writeEdgeVtuFiles(apf::Mesh* m, int n, const char* prefix)
{
  std::stringstream ss;
  ss << prefix << PCU_Comm_Self() << "_"
     << m->getShape()->getOrder()
     << getSuffix(apf::Mesh::EDGE) << ".vtu";
  int count = countOwnedEntitiesOfType(m,apf::Mesh::EDGE);
  int nCells = count*n;
  int nPoints = count*(n+1);

  std::string fileName = ss.str();
  std::ofstream file(fileName.c_str());
  assert(file.is_open());

  writeStart(file,nPoints,nCells);

  file << "<Points>\n";
  file << "<DataArray type=\"Float64\" Name=\"coordinates\" "
      "NumberOfComponents=\"3\" format=\"ascii\">\n";

  apf::MeshIterator* it = m->begin(1);
  apf::MeshEntity* e;
  apf::Vector3 p,pt;
  while ((e = m->iterate(it))) {
    if(!m->isOwned(e)) continue;
    apf::Element* elem =
        apf::createElement(m->getCoordinateField(),e);
    for (int i = 0; i <= n; ++i){
      p[0] = 2.*i/n-1.;
      apf::getVector(elem,p,pt);
      writePoint(file,pt);
    }
    apf::destroyElement(elem);
  }
  m->end(it);
  file << "</DataArray>\n";
  file << "</Points>\n";
  writeCells(file,apf::Mesh::EDGE,count,n,nCells);
  writeEnd(file);
}

static void writeTriangleVtuFiles(apf::Mesh* m, int n, const char* prefix)
{
  std::stringstream ss;
  ss << prefix << PCU_Comm_Self() << "_"
     << m->getShape()->getOrder()
     << getSuffix(apf::Mesh::TRIANGLE) << ".vtu";

  int count = countOwnedEntitiesOfType(m,apf::Mesh::TRIANGLE);
  int nPoints = count*(n+1)*(n+2)/2;
  int nCells = count*n*n;

  std::string fileName = ss.str();
  std::ofstream file(fileName.c_str());
  assert(file.is_open());

  writeStart(file,nPoints,nCells);

  file << "<Points>\n";
  file << "<DataArray type=\"Float64\" Name=\"coordinates\" "
      "NumberOfComponents=\"3\" format=\"ascii\">\n";

  apf::MeshIterator* it = m->begin(2);
  apf::MeshEntity* e;
  apf::Vector3 p,pt;
  while ((e = m->iterate(it))) {
    if(!m->isOwned(e)) continue;
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
  }
  m->end(it);
  file << "</DataArray>\n";
  file << "</Points>\n";
  writeCells(file,apf::Mesh::TRIANGLE,count,n,nCells);
  if(m->getShape()->getOrder() > 1)
    writeTriJacobianData(file,m,n);
  writeEnd(file);
}

static void writeTetVtuFiles(apf::Mesh* m, int n, const char* prefix)
{
  std::stringstream ss;
  ss << prefix << PCU_Comm_Self() << "_"
     << m->getShape()->getOrder()
     << getSuffix(apf::Mesh::TET) << ".vtu";
  std::string fileName = ss.str();
  std::ofstream file(fileName.c_str());
  assert(file.is_open());

  int count = countOwnedEntitiesOfType(m,apf::Mesh::TET);
  int nPoints = count*4*(n+1)*(n+1)*(n+1);
  int nCells = count*4*n*n*n;

  // first initializing with end points
  apf::Vector3 params[15] = {apf::Vector3(0,0,0),apf::Vector3(1,0,0),
      apf::Vector3(0,1,0),apf::Vector3(0,0,1),apf::Vector3(0,0,0),
      apf::Vector3(0,0,0),apf::Vector3(0,0,0),apf::Vector3(0,0,0),
      apf::Vector3(0,0,0),apf::Vector3(0,0,0),apf::Vector3(0,0,0),
      apf::Vector3(0,0,0),apf::Vector3(0,0,0),apf::Vector3(0,0,0),
      apf::Vector3(0.25,0.25,0.25)};

  for(int i = 0; i < 6; ++i)
    params[i+4] = params[apf::tet_edge_verts[i][0]]*0.5
      + params[apf::tet_edge_verts[i][1]]*0.5;
  for(int i = 0; i < 4; ++i)
    params[i+10] = params[apf::tet_tri_verts[i][0]]*1./3.
      + params[apf::tet_tri_verts[i][1]]*1./3.
      + params[apf::tet_tri_verts[i][2]]*1./3.;

  int hex[4][8] = {{0,4,10,6,7,11,14,13},{1,5,10,4,8,12,14,11},
      {2,6,10,5,9,13,14,12},{9,13,14,12,3,7,11,8}};

  writeStart(file,nPoints,nCells);

  file << "<Points>\n";
  file << "<DataArray type=\"Float64\" Name=\"coordinates\" "
      "NumberOfComponents=\"3\" format=\"ascii\">\n";

  apf::MeshIterator* it = m->begin(3);
  apf::MeshEntity* e;
  apf::Vector3 xi,p,pt;
  apf::NewArray<double> values;
  apf::EntityShape* shape = apf::getLagrange(1)->getEntityShape(apf::Mesh::HEX);

  while ((e = m->iterate(it))) {
    if(!m->isOwned(e)) continue;
    apf::Element* elem =
        apf::createElement(m->getCoordinateField(),e);
    for(int h = 0; h < 4; ++h){
      for (int k = 0; k <= n; ++k){
        xi[2] = 2.*k/n - 1.;
        for (int j = 0; j <= n; ++j){
          xi[1] = 2.*j/n - 1.;
          for (int i = 0; i <= n; ++i){
            xi[0] = 2.*i/n - 1.;
            shape->getValues(0, 0, xi, values);
            p.zero();
            for(int l = 0; l < 8; ++l)
              p += params[hex[h][l]]*values[l];
            apf::getVector(elem,p,pt);
            writePoint(file,pt);
          }
        }
      }
    }
    apf::destroyElement(elem);
  }

  file << "</DataArray>\n";
  file << "</Points>\n";
  writeCells(file,apf::Mesh::TET,count,n,nCells);
  writeTetJacobianData(file,m,n);
  writeEnd(file);
  m->end(it);
}

void writeControlPointVtuFiles(apf::Mesh* m, const char* prefix)
{
  if (!PCU_Comm_Self())
    writePvtuFile(prefix,m,apf::Mesh::VERTEX);

  std::stringstream ss;
  ss << prefix << PCU_Comm_Self() << "_"
     << m->getShape()->getOrder()
     << getSuffix(apf::Mesh::VERTEX) << ".vtu";

  int nPoints = 0;
    for (int t = 0; t < apf::Mesh::TYPES; ++t)
      nPoints += m->getShape()->countNodesOn(t)
      *countOwnedEntitiesOfType(m,t);

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
      if(!m->isOwned(e)) continue;
      for(int i = 0; i < m->getShape()->countNodesOn(t); ++i){
        m->getPoint(e,i,pt);
        writePoint(file,pt);
      }
    }
    m->end(it);
  }
  file << "</DataArray>\n";
  file << "</Points>\n";
  writeCells(file,apf::Mesh::VERTEX,nPoints,nPoints,nPoints);
  writeEnd(file);
}

void writeCurvedVtuFiles(apf::Mesh* m, int type, int n, const char* prefix)
{
  if (!PCU_Comm_Self())
    writePvtuFile(prefix,m,type);

  PCU_Barrier();
  switch (type) {
    case apf::Mesh::EDGE:
      writeEdgeVtuFiles(m,n,prefix);
      break;
    case apf::Mesh::TRIANGLE:
      writeTriangleVtuFiles(m,n,prefix);
      break;
    case apf::Mesh::TET:
      writeTetVtuFiles(m,n/2+1,prefix);
      break;
    default:
      break;
  }
  PCU_Barrier();
}

} //namespace crv
