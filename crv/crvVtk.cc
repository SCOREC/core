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

static int countBoundaryEdges(apf::Mesh* m)
{
  int n = 0;
  apf::MeshIterator* it = m->begin(1);
  apf::MeshEntity* e;
  while ((e = m->iterate(it))) {
    if (m->getModelType(m->toModel(e)) != m->getDimension())
      n++;
  }
  m->end(it);
  return n;
}

static int countBoundaryFaces(apf::Mesh* m, bool countEdges)
{
  int n = 0;
  apf::MeshIterator* it = m->begin(2);
  apf::MeshEntity* e;
  apf::MeshEntity* edges[3];
  while ((e = m->iterate(it))) {
    m->getDownward(e,1,edges);
    if((countEdges &&
      (m->getModelType(m->toModel(edges[0])) != m->getDimension() ||
       m->getModelType(m->toModel(edges[1])) != m->getDimension() ||
       m->getModelType(m->toModel(edges[2])) != m->getDimension())) ||
       m->getModelType(m->toModel(e)) != m->getDimension())
      n++;
  }
  m->end(it);
  return n;
}

static void writeEdgeConnectivity(std::ostream& file, apf::Mesh* m, int nSplits)
{
  file << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
  apf::MeshEntity* e;
  apf::MeshIterator* it = m->begin(1);
  int num = 0;
  while ((e = m->iterate(it)))
  {
    if (m->getModelType(m->toModel(e)) != m->getDimension()) {
      for (int i=0; i < nSplits; ++i)
        file << num+i << ' ' << num+i+1 << '\n';
      num += nSplits+1;
    } else {
      file << num << ' ' << num+1 << '\n';
      num += 2;
    }
  }
  m->end(it);
  file << "</DataArray>\n";
}
static void writeFaceConnectivity(std::ostream& file, apf::Mesh* m, int n)
{
  file << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
  apf::MeshEntity* e;
  apf::MeshEntity* edges[3];
  apf::MeshIterator* it = m->begin(2);
  int num = 0;
  while ((e = m->iterate(it)))
  {
    m->getDownward(e,1,edges);
    if(m->getModelType(m->toModel(edges[0])) != m->getDimension() ||
       m->getModelType(m->toModel(edges[1])) != m->getDimension() ||
       m->getModelType(m->toModel(edges[2])) != m->getDimension() ||
       m->getModelType(m->toModel(e)) != m->getDimension()) {
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

static void writeOffsets(std::ostream& file, int d, int nCells)
{
  file << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
  int o = 0;
  for (int i=0; i < nCells; ++i){
    o += d+1;
    file << o << '\n';
  }
  file << "</DataArray>\n";
}

static void writeTypes(std::ostream& file, int d, int nCells)
{
  file << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
  static int vtkTypes[4] = {1,3,5,10};
  for (int i=0; i < nCells; ++i)
    file << vtkTypes[d] << '\n';
  file << "</DataArray>\n";
}

static void writePoint(std::ostream& file, apf::Vector3 & pt)
{
  for (int j=0; j < 3; ++j)
    file << pt[j] << ' ';
  file << '\n';
}

/*
This will subdivide curved edges and faces and output two files, one with each
by splitting an edge into nSplit edges and nSplit*nSplit faces

The points for edges are straightforward, for faces they are connected from the
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
...
4 1 0
followed by the last 3 (n*(n-1)/2)
5 8 7
2 6 5
1 5 4
*/

static void writeEdgeVtuFiles(apf::Mesh* m, int nSplit, const char* prefix)
{
  std::stringstream ss;
  ss << prefix << PCU_Comm_Self() << "_"
     << m->getShape()->getOrder()
      << "_edges.vtu";

  int nBoundaryEnts = countBoundaryEdges(m);

  int nPoints = 2*m->count(1) + nBoundaryEnts*(nSplit-1);
  int nCells = m->count(1) + nBoundaryEnts*(nSplit-1);

  std::string fileName = ss.str();
  std::ofstream file(fileName.c_str());
  assert(file.is_open());
  file << "<VTKFile type=\"UnstructuredGrid\">\n";
  file << "<UnstructuredGrid>\n";
  file << "<Piece NumberOfPoints=\"" << nPoints;
  file << "\" NumberOfCells=\"" << nCells;
  file << "\">\n";
  file << "<Points>\n";
  file << "<DataArray type=\"Float64\" Name=\"coordinates\" "
      "NumberOfComponents=\"3\" format=\"ascii\">\n";

  apf::MeshIterator* it = m->begin(1);
  apf::MeshEntity* e;
  apf::MeshEntity* v[2];
  apf::Vector3 pa,pt;
  while ((e = m->iterate(it))) {
    if (m->getModelType(m->toModel(e)) != m->getDimension()){
      apf::Element* elem =
          apf::createElement(m->getCoordinateField(),e);
      for (int i = 0; i <= nSplit; ++i){
        pa[0] = 2.*i/nSplit-1.;
        apf::getVector(elem,pa,pt);
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
  file << "<Cells>\n";

  writeEdgeConnectivity(file,m,nSplit);
  writeOffsets(file,1,nCells);
  writeTypes(file,1,nCells);

  file << "</Cells>\n";
  file << "</Piece>\n";
  file << "</UnstructuredGrid>\n";
  file << "</VTKFile>\n";
}

static void writeFaceVtuFiles(apf::Mesh* m, int nSplit, const char* prefix)
{
  std::stringstream ss;
  ss << prefix << PCU_Comm_Self() << "_"
     << m->getShape()->getOrder()
     << "_faces.vtu";

  int nBoundaryEnts = countBoundaryFaces(m,true);

  int nPoints = 3*m->count(2) + nBoundaryEnts*((nSplit+1)*(nSplit+2)/2-3);
  int nCells = m->count(2) + nBoundaryEnts*(nSplit*nSplit-1);

  std::string fileName = ss.str();
  std::ofstream file(fileName.c_str());
  assert(file.is_open());
  file << "<VTKFile type=\"UnstructuredGrid\">\n";
  file << "<UnstructuredGrid>\n";
  file << "<Piece NumberOfPoints=\"" << nPoints;
  file << "\" NumberOfCells=\"" << nCells;
  file << "\">\n";
  file << "<Points>\n";
  file << "<DataArray type=\"Float64\" Name=\"coordinates\" "
      "NumberOfComponents=\"3\" format=\"ascii\">\n";

  apf::MeshIterator* it = m->begin(2);
  apf::MeshEntity* e;
  apf::MeshEntity* down[3];
  apf::Vector3 pa,pt;
  while ((e = m->iterate(it))) {

    m->getDownward(e,1,down);
    if(m->getModelType(m->toModel(down[0])) != m->getDimension() ||
       m->getModelType(m->toModel(down[1])) != m->getDimension() ||
       m->getModelType(m->toModel(down[2])) != m->getDimension() ||
       m->getModelType(m->toModel(e)) != m->getDimension()){
      apf::Element* elem =
          apf::createElement(m->getCoordinateField(),e);
      for (int j = 0; j <= nSplit; ++j){
        pa[1] = 1.*j/nSplit;
        for (int i = 0; i <= nSplit-j; ++i){
          pa[0] = 1.*i/nSplit;
          apf::getVector(elem,pa,pt);
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
  file << "<Cells>\n";

  writeFaceConnectivity(file,m,nSplit);
  writeOffsets(file,2,nCells);
  writeTypes(file,2,nCells);

  file << "</Cells>\n";
  file << "</Piece>\n";
  file << "</UnstructuredGrid>\n";
  file << "</VTKFile>\n";
}
static void writeControlPointVtuFiles(apf::Mesh* m, const char* prefix)
{
  std::stringstream ss;
  ss << prefix << PCU_Comm_Self() << "_"
     << m->getShape()->getOrder()
     << "_controlPoints.vtu";

  int type = m->getDimension() == 2 ? apf::Mesh::EDGE : apf::Mesh::TRIANGLE;
  int nBase = m->getDimension() == 2 ? 2 : 3;
  int nNodes = m->getShape()->getEntityShape(type)->countNodes();

  int nBoundaryEnts = m->getDimension() == 2 ? countBoundaryEdges(m) :
      countBoundaryFaces(m,false);
  int nPoints = nBase*m->count(2) + nBoundaryEnts*(nNodes-nBase);

  std::string fileName = ss.str();
  std::ofstream file(fileName.c_str());
  assert(file.is_open());
  file << "<VTKFile type=\"UnstructuredGrid\">\n";
  file << "<UnstructuredGrid>\n";
  file << "<Piece NumberOfPoints=\"" << nPoints;
  file << "\" NumberOfCells=\"" << nPoints;
  file << "\">\n";
  file << "<Points>\n";
  file << "<DataArray type=\"Float64\" Name=\"coordinates\" "
      "NumberOfComponents=\"3\" format=\"ascii\">\n";
  apf::MeshIterator* it = m->begin(m->getDimension()-1);
  apf::MeshEntity* e;
  apf::MeshEntity* v[3];
  apf::Vector3 pa,pt;
  apf::NewArray<apf::Vector3> nodes;

  while ((e = m->iterate(it))) {
    if(m->getModelType(m->toModel(e)) != m->getDimension()){
      apf::Element* elem =
          apf::createElement(m->getCoordinateField(),e);
      apf::getVectorNodes(elem,nodes);
      for(int i = 0; i < nNodes; ++i)
        writePoint(file,nodes[i]);

      apf::destroyElement(elem);

    } else {
      int nVerts = m->getDownward(e,0,v);
      for(int i = 0; i < nVerts; ++i)
        writePoint(file,pt);

    }
  }
  m->end(it);
  file << "</DataArray>\n";
  file << "</Points>\n";
  file << "<Cells>\n";

  file << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
  for(int i = 0; i < nPoints; ++i)
    file << i << '\n';
  file << "</DataArray>\n";  writeOffsets(file,0,nPoints);
  writeTypes(file,0,nPoints);

  file << "</Cells>\n";
  file << "</Piece>\n";
  file << "</UnstructuredGrid>\n";
  file << "</VTKFile>\n";
}

void writeCurvedVtuFiles(apf::Mesh* m, int n, const char* prefix)
{
  writeEdgeVtuFiles(m,n,prefix);
  writeFaceVtuFiles(m,n,prefix);
  writeControlPointVtuFiles(m,prefix);
}

} //namespace crv
