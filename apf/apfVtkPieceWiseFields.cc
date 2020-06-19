/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include <PCU.h>
#include <lionBase64.h>
#include <lionCompress.h>
#include "apfMesh.h"
#include "apfNumbering.h"
#include "apfNumberingClass.h"
#include "apfShape.h"
#include "apfFieldData.h"
#include <sstream>
#include <fstream>
#include <pcu_util.h>
#include <lionPrint.h>
#include <cstdlib>
#include <stdint.h>
#include <vector>
#include <apfVtk.h>

// === includes for safe_mkdir ===
#include <reel.h>
#include <sys/types.h> /*required for mode_t for mkdir on some systems*/
#include <sys/stat.h> /*using POSIX mkdir call for SMB "foo/" path*/
#include <errno.h> /* for checking the error from mkdir */
// ===============================

namespace apf {

static void safe_mkdir(const char* path)
{
  mode_t const mode = S_IRWXU|S_IRGRP|S_IXGRP|S_IROTH|S_IXOTH;
  int err;
  errno = 0;
  err = mkdir(path, mode);
  if (err != 0 && errno != EEXIST)
  {
    reel_fail("MDS: could not create directory \"%s\"\n", path);
  }
}

/* Paraview/VTK has trouble with sub-normal double precision floating point
 * ASCII values.
 *
 * http://www.paraview.org/Bug/view.php?id=15925
 *
 * This function exists to cast "double" to "float"
 * before writing it to file, the others are to maintain the
 * templated design of writeCornerCoords and others */


static std::string getPieceFileName(int id)
{
  std::stringstream ss;
  ss << id << ".vtk";
  return ss.str();
}

static std::string getFileNameAndPathVtu(const char* prefix,
    std::string fileName,
    int id)
{
  std::stringstream ss;
  ss << prefix << '/' << prefix << '_' << id << '_' << fileName;
  return ss.str();
}

static void writeNedelecVtkFile(const char* prefix, Mesh* m,
    std::vector<std::string> writeFields)
{
  double t0 = PCU_Time();

  // get the number of points on this part
  MeshEntity* e;
  MeshIterator* it = m->begin(m->getDimension());
  int point_count = 0;
  int cell_count = 0;
  bool isSimplexMesh = true;
  while( (e = m->iterate(it)) ) {
    int type = m->getType(e);
    if (!apf::isSimplex(type)) {
      isSimplexMesh = false;
      break;
    }
    point_count += Mesh::adjacentCount[type][0];
    cell_count++;
  }
  m->end(it);

  PCU_ALWAYS_ASSERT_VERBOSE(isSimplexMesh,
      "writeNedelecVtk only implemented for all simplex meshes");

  static Vector3 xis[4] = {
    apf::Vector3(0., 0., 0.),
    apf::Vector3(1., 0., 0.),
    apf::Vector3(0., 1., 0.),
    apf::Vector3(0., 0., 1.)
  };

  std::string fileName = getPieceFileName(PCU_Comm_Self());
  std::string fileNameAndPath = getFileNameAndPathVtu(prefix, fileName, PCU_Comm_Self());
  std::stringstream buf;
  buf <<
    "# vtk DataFile Version 3.0\n"
    "Made by SCOREC/core\n"
    "ASCII\n"
    "DATASET UNSTRUCTURED_GRID\n";
  buf << "POINTS " << point_count << " double\n";

  it = m->begin(m->getDimension());
  while( (e = m->iterate(it)) ) {
    int type = m->getType(e);
    int n = Mesh::adjacentCount[type][0];
    MeshElement* me = createMeshElement(m, e);
    for(int i=0; i<n; i++) {
      Vector3 p;
      mapLocalToGlobal(me, xis[i], p);
      buf << p[0] << ' ' << p[1] << ' ' << p[2] << '\n';
    }
    destroyMeshElement(me);
  }
  m->end(it);

  buf << "CELLS " << cell_count << ' ' << point_count + cell_count << '\n';

  it = m->begin(m->getDimension());
  int N = 0;
  while( (e = m->iterate(it)) ) {
    int type = m->getType(e);
    int n = Mesh::adjacentCount[type][0];
    buf << n;
    for(int i=0; i<n; i++) {
      buf << ' ' << N;
      N++;
    }
    buf << '\n';
  }
  m->end(it);

  buf << "CELL_TYPES " << cell_count << '\n';
  it = m->begin(m->getDimension());
  while( (e = m->iterate(it)) ) {
    int type = m->getType(e);
    int vtk_type = -1;
    switch (type) {
      case Mesh::TRIANGLE: vtk_type = 5; break;
      case Mesh::TET:      vtk_type = 10; break;
      default:
        PCU_ALWAYS_ASSERT_VERBOSE(0,
            "only TRIANGLE and TET supported");
        break;
    }
    buf << vtk_type << '\n';
  }
  m->end(it);

  buf << "CELL_DATA " << cell_count << '\n';
  buf << "SCALARS part_id int\n";
  buf << "LOOKUP_TABLE default\n";
  it = m->begin(m->getDimension());
  while( (e = m->iterate(it)) ) {
    buf << PCU_Comm_Self() << '\n';
  }
  m->end(it);

  buf << "POINT_DATA " << point_count << '\n' << std::flush;

  for(std::size_t i = 0; i < writeFields.size(); i++) {
    Field* f = m->findField(writeFields[i].c_str());
    const char* fName = f->getName();
    buf << "VECTORS " << fName << " double\n";

    it = m->begin(m->getDimension());
    while( (e = m->iterate(it)) ) {
      int type = m->getType(e);
      int n = Mesh::adjacentCount[type][0];
      MeshElement* me = createMeshElement(m, e);
      Element* el = createElement(f, me);
      for(int j=0; j<n; j++) {
	Vector3 fVal;
	getVector(el, xis[j], fVal);
	buf << fVal[0] << ' ' << fVal[1] << ' ' << fVal[2] << '\n';
      }
      destroyElement(el);
      destroyMeshElement(me);
    }
    m->end(it);
  }
  double t1 = PCU_Time();
  if (!PCU_Comm_Self())
  {
    lion_oprint(1,"writeVtuFile into buffers: %f seconds\n", t1 - t0);
  }
  { //block forces std::ofstream destructor call
    std::ofstream file(fileNameAndPath.c_str());
    PCU_ALWAYS_ASSERT(file.is_open());
    file << buf.rdbuf();
  }
  double t2 = PCU_Time();
  if (!PCU_Comm_Self())
  {
    lion_oprint(1,"writeNedelecVtkFile buffers to disk: %f seconds\n", t2 - t1);
  }
}

std::vector<std::string> populateNedelecFields(Mesh* m)
{
  std::vector<std::string> writeFields;
  for (int i=0; i < m->countFields(); ++i)
  {
    Field* f = m->getField(i);
    std::string fsName = f->getShape()->getName();
    if (fsName == std::string("Nedelec"))
      writeFields.push_back(f->getName());
  }
  return writeFields;
}

void writeNedelecVtkFiles(const char* prefix, Mesh* m)
{
  // bool isWritingBinary = true;
  std::vector<std::string> writeFields = populateNedelecFields(m);
  if (writeFields.size() == 0) return;
  safe_mkdir(prefix);
  writeNedelecVtkFile(prefix, m, writeFields);
}

}
