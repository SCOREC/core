#include "apfMPAS.h"
#include <apfMesh2.h>
#include <apfNumbering.h>
#include <diffMC/src/parma_ghostOwner.h>
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
  return 0;
}

void numberInitialNodes(MpasFile& in, apf::Mesh2* out,
    apf::MeshEntity** v)
{
  apf::Numbering* n = apf::createNumbering(out, "mpas_id", out->getShape(), 1);
  for (int i = 0; i < in.nodes; ++i)
    apf::number(n, v[i], 0, 0, i);
}

void buildInitialMesh(MpasFile& in, apf::Mesh2* out)
{
  apf::MeshEntity** v = new apf::MeshEntity* [in.nodes];
/* fake classification for now, deriveMdsModel will fix this up */
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
/* zero indices (-1 after conversion) indicate that an MPAS
   vertex does not have an adjacent cell in that direction
   (it is a boundary vertex). In our dual mesh, we simply
   omit that dual element. This is correct and should not
   cause problems. */
      if (vi < 0)
        goto bad_vi;
      assert(vi < in.nodes);
      ev[j] = v[vi];
    }
    apf::buildElement(out, c, t, ev);
bad_vi:
    continue;
  }
  numberInitialNodes(in, out, v);
  delete [] v;
}

void removeIsolatedNodes(apf::Mesh2* m)
{
/* another problem with the dual: MPAS cells without adjacent
   cells become vertices with no adjacent elements in the
   dual mesh, which our codes don't handle.
   We have to give up on these cells, but tell the user about it.
   We also need to be aware of these omitted cells going forward */
  //apf::Numbering* nums = apf::createNumbering(m, "mpas_id", m->getShape(), 1);
  apf::MeshIterator* it = m->begin(0);
  int n = 0;
  apf::MeshEntity* e;
  while ((e = m->iterate(it)))
    if ( ! m->countUpward(e)) {
      /*int num =getNumber(nums, e, 0, 0);
        fprintf(stdout, "Missing vertex with number %d\n",num);
      */
      m->destroy(e);
      ++n;
    }
  m->end(it);
  if (n)
    fprintf(stderr,
        "warning: removed %d isolated nodes from MPAS dual mesh\n",
        n);
}

void loadMpasMesh(apf::Mesh2* m, const char* filename)
{
  assert(PCU_Comm_Peers() == 1);
  MpasFile f;
  readMpasFile(f, filename);
  buildInitialMesh(f, m);
  removeIsolatedNodes(m);
}


void writeMpasAssignments(apf::Mesh2* m, const char* filename) {
  NcFile in(filename, NcFile::read);
  /* this is the dual of a hexagonal mesh,
     hence the reversing of terms */
  int numMpasVtx = readDim(in, "nCells");
  
  apf::Numbering* n = apf::createNumbering(m, "mpas_id", m->getShape(), 1);

  /* collect N/#parts contiguous vertex assignments on each process
     (and deal with any remainders) */
  int numPerPart = numMpasVtx / PCU_Comm_Peers();
  apf::MeshIterator* itr = m->begin(0);
  apf::MeshEntity* e;
  int size;
  if (PCU_Comm_Self()==PCU_Comm_Peers()-1)
    size=numMpasVtx%numPerPart+numPerPart;
  else
    size=numPerPart;
  std::vector<int> vtxs(size,-1);
  PCU_Comm_Begin();
  while ((e = m->iterate(itr))) {
    if (!parma::isOwned(m, e))
      continue;
    int num = getNumber(n, e, 0, 0);
    int target = num / numPerPart;
    if (target>=PCU_Comm_Peers())
      target = PCU_Comm_Peers()-1;
    if (target == PCU_Comm_Self())
      vtxs[num%size] = PCU_Comm_Self();
    else
      PCU_COMM_PACK(target,num);
  }
  m->end(itr);

  PCU_Comm_Send();
  while (PCU_Comm_Receive()) {
    int owner =PCU_Comm_Sender();
    int num;
    PCU_COMM_UNPACK(num);
    vtxs[num%size]=owner;
  }

  // assign missing vertices to a random part id
  int count = 0;
  for (int i = 0;
       (i < size);
       i++)
    if (vtxs[i]==-1) {
      vtxs[i] = 0; //to be random
      //fprintf(stdout,"missing vertex %d\n",i+numPerPart*PCU_Comm_Self());
      count++;
    }
  fprintf(stdout,"missing vertices found %d\n",count);
  
  // use MPI IO to write the contiguous blocks to a single graph.info.part.<#parts> file
  // see https://gist.github.com/cwsmith/166d5beb400f3a8136f7 and the comments
  double startTime=MPI_Wtime();
  FILE* file;
  char name[32];
  
  sprintf(name,"graph.info.part.%d",PCU_Comm_Peers());
  file = fopen(name, "w");
  fseek(file,numPerPart*PCU_Comm_Self()*16,SEEK_SET);
  for (int i=0;i<size;i++)
    fprintf(file,"%-15d\n",vtxs[i]);
  
  fclose(file);
  double totalTime= MPI_Wtime()-startTime; 
  PCU_Max_Doubles(&totalTime,1);
  if (!PCU_Comm_Self())
    fprintf(stdout,"File writing time: %f seconds\n", totalTime); 
}

}
