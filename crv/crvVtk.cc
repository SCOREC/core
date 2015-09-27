/*
 * Copyright 2015 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "crv.h"
#include "PCU.h"
#include "apfDynamicVector.h"
#include "apfFieldData.h"
#include <sstream>
#include <fstream>
#include <cassert>

namespace crv {

class HasAll : public apf::FieldOp
{
  public:
      virtual bool inEntity(apf::MeshEntity* e)
      {
        if (!f->getData()->hasEntity(e) && f->countNodesOn(e))
          ok = false;
        return false;
      }
      bool run(apf::FieldBase* f_)
      {
        f = f_;
        ok = true;
        this->apply(f);
        return ok;
      }
  private:
    bool ok;
    apf::FieldBase* f;
};

static bool isPrintable(apf::FieldBase* f)
{
  HasAll op;
  return op.run(f);
}

static void describeArray(
    std::ostream& file,
    const char* name,
    int type,
    int size)
{
  file << "type=\"";
  const char* typeNames[3] = {"Float64","Int32","Int64"};
  file << typeNames[type];
  file << "\" Name=\"" << name;
  file << "\" NumberOfComponents=\"" << size;
  file << "\" format=\"ascii\"";
}

static void writeDataHeader(std::ostream& file, const char* name,
    int type, int size)
{
  file << "<DataArray ";
  describeArray(file,name,type,size);
  file << ">\n";
}

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

static void writeEdgeJacobianDet(std::ostream& file, apf::Mesh* m, int n)
{
  file << "<DataArray type=\"Float64\" Name=\"detJacobian\" "
       << "NumberOfComponents=\"1\" format=\"ascii\">\n";

  apf::MeshIterator* it = m->begin(1);
  apf::MeshEntity* e;
  apf::Vector3 p,pt;

  apf::Matrix3x3 J;

  while ((e = m->iterate(it))) {
    if(!m->isOwned(e)) continue;
    apf::MeshElement* me = apf::createMeshElement(m,e);
    for (int i = 0; i <= n; ++i){
      p[0] = 2.*i/n-1.;
      apf::getJacobian(me,p,J);
      double detJ = apf::getJacobianDeterminant(J,1);
      file << detJ << '\n';
    }
    apf::destroyMeshElement(me);
  }
  m->end(it);

  file << "</DataArray>\n";
}

static void writeTriJacobianDet(std::ostream& file, apf::Mesh* m, int n)
{
  file << "<DataArray type=\"Float64\" Name=\"detJacobian\" "
       << "NumberOfComponents=\"1\" format=\"ascii\">\n";

  apf::MeshIterator* it = m->begin(2);
  apf::MeshEntity* e;
  apf::Vector3 p;

  apf::Matrix3x3 J;
  double detJ;
  bool isValid = true;
  while ((e = m->iterate(it))) {
    if(!m->isOwned(e) || m->getType(e) != apf::Mesh::TRIANGLE) continue;
    apf::MeshElement* me = apf::createMeshElement(m,e);
    for (int j = 0; j <= n; ++j){
      p[1] = 1.*j/n;
      for (int i = 0; i <= n-j; ++i){
        p[0] = 1.*i/n;
        apf::getJacobian(me,p,J);
        if(m->getDimension() == 3){
          detJ = apf::getJacobianDeterminant(J,2);
        } else {
          detJ = J[0][0]*J[1][1]-J[1][0]*J[0][1];
        }
        file << detJ << '\n';
        if(isValid && detJ < 0.){
          apf::Vector3 pt;
          apf::getVector((apf::Element*)me,p,pt);
          std::stringstream ss;
          ss << "warning: Tri Jacobian Determinant is negative,  " << detJ
             << '\n';
          std::string s = ss.str();
          fprintf(stderr, "%s", s.c_str());
          isValid = false;
        }
      }
    }
    apf::destroyMeshElement(me);
  }
  m->end(it);

  file << "</DataArray>\n";
}

static void writeTetJacobianDet(std::ostream& file, apf::Mesh* m, int n)
{
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

  bool isValid = true;
  while ((e = m->iterate(it))) {
    if(!m->isOwned(e) || m->getType(e) != apf::Mesh::TET) continue;
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
            if(isValid && detJ < 0.){
              apf::Vector3 pt;
              apf::getVector((apf::Element*)me,p,pt);
              std::stringstream ss;
              ss << "warning: Tet Jacobian Determinant is negative,  " << detJ
                 << '\n';
              std::string s = ss.str();
              fprintf(stderr, "%s", s.c_str());
              isValid = false;
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
}

static void writeJacobianDet(std::ostream& file, apf::Mesh* m, int type, int n)
{
  switch (type) {
    case apf::Mesh::EDGE:
      writeEdgeJacobianDet(file,m,n);
      break;
    case apf::Mesh::TRIANGLE:
      writeTriJacobianDet(file,m,n);
      break;
    case apf::Mesh::TET:
      writeTetJacobianDet(file,m,n);
      break;
    default:
      break;
  }
}

static void writeEdgeNodalField(std::ostream& file, int n, apf::Field* f)
{
  int nc = f->countComponents();
  writeDataHeader(file,f->getName(),f->getScalarType(),nc);
  apf::Mesh* m = f->getMesh();
  apf::MeshIterator* it = m->begin(1);
  apf::MeshEntity* e;
  apf::Vector3 p;
  apf::DynamicVector v(nc);
  while ((e = m->iterate(it))) {
    if(!m->isOwned(e)) continue;
    apf::Element* elem = apf::createElement(f,e);
    for (int i = 0; i <= n; ++i){
      p[0] = 2.*i/n-1.;
      apf::getComponents(elem,p,&v[0]);
      for (int j = 0; j < nc; ++j)
        file << v[j] << ' ';
      file << '\n';
    }
    apf::destroyElement(elem);
  }
  m->end(it);
  file << "</DataArray>\n";
}

static void writeTriNodalField(std::ostream& file, int n, apf::Field* f)
{
  int nc = f->countComponents();
  writeDataHeader(file,f->getName(),f->getScalarType(),nc);
  apf::Mesh* m = f->getMesh();
  apf::MeshIterator* it = m->begin(2);
  apf::MeshEntity* e;
  apf::Vector3 p;
  apf::DynamicVector v(nc);
  while ((e = m->iterate(it))) {
    if(!m->isOwned(e) || m->getType(e) != apf::Mesh::TRIANGLE) continue;
    apf::Element* elem =
        apf::createElement(f,e);
    for (int j = 0; j <= n; ++j){
      p[1] = 1.*j/n;
      for (int i = 0; i <= n-j; ++i){
        p[0] = 1.*i/n;
        apf::getComponents(elem,p,&v[0]);
        for (int j = 0; j < nc; ++j)
          file << v[j] << ' ';
        file << '\n';
      }
    }
    apf::destroyElement(elem);
  }
  m->end(it);
  file << "</DataArray>\n";
}

static void writeTetNodalField(std::ostream& file, int n, apf::Field* f)
{
  int nc = f->countComponents();
  writeDataHeader(file,f->getName(),f->getScalarType(),nc);
  apf::Mesh* m = f->getMesh();
  apf::MeshIterator* it = m->begin(3);
  apf::MeshEntity* e;
  apf::Vector3 xi,p;
  apf::DynamicVector v(nc);
  apf::NewArray<double> values;
  apf::EntityShape* shape = apf::getLagrange(1)->getEntityShape(apf::Mesh::HEX);

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

  while ((e = m->iterate(it))) {
    if(!m->isOwned(e) || m->getType(e) != apf::Mesh::TET) continue;
    apf::Element* elem =
        apf::createElement(f,e);
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
            apf::getComponents(elem,p,&v[0]);
            for (int j = 0; j < nc; ++j)
              file << v[j] << ' ';
            file << '\n';
          }
        }
      }
    }
    apf::destroyElement(elem);
  }
  m->end(it);

  file << "</DataArray>\n";
}

static void writeNodalField(std::ostream& file, int type, int n,
    apf::Field* f)
{
  switch (type) {
    case apf::Mesh::EDGE:
      writeEdgeNodalField(file,n,f);
      break;
    case apf::Mesh::TRIANGLE:
      writeTriNodalField(file,n,f);
      break;
    case apf::Mesh::TET:
      writeTetNodalField(file,n,f);
      break;
    default:
      break;
  }
}

static void writePDataArray(
    std::ostream& file,
    apf::FieldBase* f)
{
  file << "<PDataArray ";
  describeArray(file,
      f->getName(),
      f->getScalarType(),
      f->countComponents());
  file << "/>\n";
}

static void writePPointData(std::ostream& file, apf::Mesh* m)
{
  for (int i=0; i < m->countFields(); ++i)
  {
    apf::Field* f = m->getField(i);
    if(isPrintable(f))
      writePDataArray(file,f);
  }
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
  writePDataArray(file,m->getCoordinateField());
  file << "</PPoints>\n";
  file << "<PPointData>\n";
  if(type == apf::Mesh::VERTEX){
    file << "<PDataArray type=\"UInt8\" Name=\"entityType\" "
         << "NumberOfComponents=\"1\" format=\"ascii\"/>\n";
  } else {
    file << "<PDataArray type=\"Float64\" Name=\"detJacobian\" "
         << "NumberOfComponents=\"1\" format=\"ascii\"/>\n";
  }
  writePPointData(file,m);
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

static void writePointData(std::ostream& file, apf::Mesh* m,
    int type, int n)
{
  for (int i=0; i < m->countFields(); ++i)
  {
    apf::Field* f = m->getField(i);
    if(isPrintable(f))
      writeNodalField(file,type,n,f);
  }
}

void writeControlPointVtuFiles(apf::Mesh* m, const char* prefix)
{
  if (!PCU_Comm_Self())
    writePvtuFile(prefix,m,apf::Mesh::VERTEX);

  PCU_Barrier();

  std::stringstream ss;
  ss << prefix << PCU_Comm_Self() << "_"
     << m->getShape()->getOrder()
     << "_interPts" << ".vtu";

  int nPoints = 0;
    for (int t = 0; t < apf::Mesh::TYPES; ++t)
      nPoints += m->getShape()->countNodesOn(t)
      *countOwnedEntitiesOfType(m,t);

  std::string fileName = ss.str();
  std::stringstream buf;

  writeStart(buf,nPoints,nPoints);
  buf << "<Points>\n";
  buf << "<DataArray type=\"Float64\" Name=\"coordinates\" "
      "NumberOfComponents=\"3\" format=\"ascii\">\n";

  for (int t = 0; t < apf::Mesh::TYPES; ++t){
    apf::MeshIterator* it = m->begin(apf::Mesh::typeDimension[t]);
    apf::MeshEntity* e;
    apf::Vector3 pt, xi;
    while ((e = m->iterate(it))) {
      if(!m->isOwned(e)) continue;
      apf::Element* elem = apf::createElement(m->getCoordinateField(),e);
      for(int i = 0; i < m->getShape()->countNodesOn(t); ++i){
        m->getShape()->getNodeXi(t,i,xi);
        apf::getVector(elem,xi,pt);
        writePoint(buf,pt);
      }
      apf::destroyElement(elem);
    }
    m->end(it);
  }
  buf << "</DataArray>\n";
  buf << "</Points>\n";
  writeCells(buf,apf::Mesh::VERTEX,nPoints,nPoints,nPoints);
  buf << "<PointData>\n";
  buf << "<DataArray type=\"UInt8\" Name=\"entityType\" "
      << "NumberOfComponents=\"1\" format=\"ascii\">\n";

  for (int t = 0; t < apf::Mesh::TYPES; ++t){
    apf::MeshIterator* it = m->begin(apf::Mesh::typeDimension[t]);
    apf::MeshEntity* e;
    apf::Vector3 pt;
    while ((e = m->iterate(it))) {
      if(!m->isOwned(e)) continue;
      for(int i = 0; i < m->getShape()->countNodesOn(t); ++i){
        buf << t << '\n';
      }
    }
    m->end(it);
  }
  buf << "</DataArray>\n";
  buf << "</PointData>\n";
  writeEnd(buf);
  {
    std::ofstream file(fileName.c_str());
    assert(file.is_open());
    file << buf.rdbuf();
  }

  PCU_Barrier();
}

void writeInterpolationPointVtuFiles(apf::Mesh* m, const char* prefix)
{
  if (!PCU_Comm_Self())
    writePvtuFile(prefix,m,apf::Mesh::VERTEX);

  PCU_Barrier();

  std::stringstream ss;
  ss << prefix << PCU_Comm_Self() << "_"
     << m->getShape()->getOrder()
     << getSuffix(apf::Mesh::VERTEX) << ".vtu";

  int nPoints = 0;
    for (int t = 0; t < apf::Mesh::TYPES; ++t)
      nPoints += m->getShape()->countNodesOn(t)
      *countOwnedEntitiesOfType(m,t);

  std::string fileName = ss.str();
  std::stringstream buf;

  writeStart(buf,nPoints,nPoints);
  buf << "<Points>\n";
  buf << "<DataArray type=\"Float64\" Name=\"coordinates\" "
      "NumberOfComponents=\"3\" format=\"ascii\">\n";

  for (int t = 0; t < apf::Mesh::TYPES; ++t){
    apf::MeshIterator* it = m->begin(apf::Mesh::typeDimension[t]);
    apf::MeshEntity* e;
    apf::Vector3 pt;
    while ((e = m->iterate(it))) {
      if(!m->isOwned(e)) continue;
      for(int i = 0; i < m->getShape()->countNodesOn(t); ++i){
        m->getPoint(e,i,pt);
        writePoint(buf,pt);
      }
    }
    m->end(it);
  }
  buf << "</DataArray>\n";
  buf << "</Points>\n";
  writeCells(buf,apf::Mesh::VERTEX,nPoints,nPoints,nPoints);
  buf << "<PointData>\n";
  buf << "<DataArray type=\"UInt8\" Name=\"entityType\" "
      << "NumberOfComponents=\"1\" format=\"ascii\">\n";

  for (int t = 0; t < apf::Mesh::TYPES; ++t){
    apf::MeshIterator* it = m->begin(apf::Mesh::typeDimension[t]);
    apf::MeshEntity* e;
    apf::Vector3 pt;
    while ((e = m->iterate(it))) {
      if(!m->isOwned(e)) continue;
      for(int i = 0; i < m->getShape()->countNodesOn(t); ++i){
        buf << t << '\n';
      }
    }
    m->end(it);
  }
  buf << "</DataArray>\n";
  buf << "</PointData>\n";
  writeEnd(buf);
  {
    std::ofstream file(fileName.c_str());
    assert(file.is_open());
    file << buf.rdbuf();
  }

  PCU_Barrier();
}

void writeCurvedVtuFiles(apf::Mesh* m, int type, int n, const char* prefix)
{
  double t0 = PCU_Time();
  if (!PCU_Comm_Self())
    writePvtuFile(prefix,m,type);

  PCU_Barrier();

  std::stringstream ss;
  ss << prefix << PCU_Comm_Self() << "_"
     << m->getShape()->getOrder()
     << getSuffix(type) << ".vtu";
  std::string fileName = ss.str();
  std::stringstream buf;

  int nPoints = 0, nCells = 0;
  int count = countOwnedEntitiesOfType(m,type);

  switch (type) {
  case apf::Mesh::EDGE:
    nCells = count*n;
    nPoints = count*(n+1);
    break;
  case apf::Mesh::TRIANGLE:
    nPoints = count*(n+1)*(n+2)/2;
    nCells = count*n*n;
    break;
  case apf::Mesh::TET:
    n = n/2+1;
    nPoints = count*4*(n+1)*(n+1)*(n+1);
    nCells = count*4*n*n*n;
    break;
  default:
    break;
  }

  writeStart(buf,nPoints,nCells);

  buf << "<Points>\n";
  writeNodalField(buf,type,n,m->getCoordinateField());
  buf << "</Points>\n";
  writeCells(buf,type,count,n,nCells);
  buf << "<PointData>\n";
  writeJacobianDet(buf,m,type,n);
  writePointData(buf,m,type,n);
  buf << "</PointData>\n";
  writeEnd(buf);

  {
    std::ofstream file(fileName.c_str());
    assert(file.is_open());
    file << buf.rdbuf();
  }

  PCU_Barrier();
  double t1 = PCU_Time();
  if (!PCU_Comm_Self())
    printf("%s vtk files %s written in %f seconds\n",
        apf::Mesh::typeName[type],prefix, t1 - t0);
}

} //namespace crv
