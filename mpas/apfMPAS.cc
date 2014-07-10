#include "apfMPAS.h"
#include <apfMesh2.h>
#include <netcdf>

namespace apf {

using namespace netCDF;

struct MpasFile
{
  ~MpasFile()
  {
    delete [] x;
    delete [] y;
    delete [] z;
    delete [] conn;
  }
  int elements;
  int nodes;
  int nodesPerElement;
  double* x;
  double* y;
  double* z;
/* conn[i * vertexDegree + j] - 1
   is the index of the j'th node
   of element i */
  int* conn;
};

static int readDim(NcFile& f, const char* name)
{
  return f.getDim(name).getSize();
}

static double* readDoubles(NcFile& f, const char* name, int n)
{
  double* a = new double[n];
  NcVar v = f.getVar(name);
  v.getVar(a);
  return a;
}

static int* readInts(NcFile& f, const char* name, int n)
{
  int* a = new int[n];
  NcVar v = f.getVar(name);
  v.getVar(a);
  return a;
}

void readMpasFile(MpasFile& out, const char* filename)
{
  NcFile in(filename, NcFile::read);
  /* this is the dual of a hexagonal mesh,
     hence the reversing of terms */
  out.nodes = readDim(in, "nCells");
  printf("%d nodes\n", out.nodes);
  out.elements = readDim(in, "nVertices");
  printf("%d elements\n", out.elements);
  out.nodesPerElement = readDim(in, "vertexDegree");
  out.x = readDoubles(in, "xCell", out.nodes);
  out.y = readDoubles(in, "yCell", out.nodes);
  out.z = readDoubles(in, "zCell", out.nodes);
  out.conn = readInts(in, "cellsOnVertex", out.elements * out.nodesPerElement);
}

static int getType(int nodesPerElement)
{
  if (nodesPerElement == 3)
    return apf::Mesh::TRIANGLE;
  if (nodesPerElement == 4)
    return apf::Mesh::QUAD;
  abort();
}

void buildMesh(MpasFile& in, apf::Mesh2* out)
{
  apf::MeshEntity** v = new apf::MeshEntity* [in.nodes];
  apf::ModelEntity* c = out->findModelEntity(2, 0);
  int t = getType(in.nodesPerElement);
  for (int i = 0; i < in.nodes; ++i)
    v[i] = out->createVert(c);
  for (int i = 0; i < in.nodes; ++i) {
    apf::Vector3 x(in.x[i], in.y[i], in.z[i]);
    out->setPoint(v[i], 0, x);
  }
  for (int i = 0; i < in.elements; ++i) {
    apf::Downward ev;
    for (int j = 0; j < in.nodesPerElement; ++j) {
      int vi = in.conn[i * in.nodesPerElement + j] - 1;
/* MPAS files seem to have connectivity numbered from 1,
   (the value nCells occurs in the indices),
   but many entries have zero indices, and whats worse
   some of them have multiple zero indices.
   this is reflected in the reader from VTK as well.
   We ditch anything with a zero (-1 after convert) index. */
      if (vi < 0)
        goto bad_vi;
      assert(vi < in.nodes);
      ev[j] = v[vi];
    }
    apf::buildElement(out, c, t, ev);
bad_vi:
    continue;
  }
  delete [] v;
}

void loadMpasMesh(apf::Mesh2* m, const char* filename)
{
  MpasFile f;
  readMpasFile(f, filename);
  buildMesh(f, m);
}

}
