#include <gmi_mesh.h>
#include <gmi_sim.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <apfConvert.h>
#include <apf.h>
#include <PCU.h>
#include <SimUtil.h>
#include <cassert>
#include <cstdlib>

#define INTERIORTAG 0
#define BOTTOMTAG 1
#define TOPTAG 2
#define PERIMETERTAG 3
#define BOTTOM_PERIMETERTAG 4
#define TOP_PERIMETERTAG 5
#define TOPFACE 1
#define BOTTOMFACE 2
#define PERIMETERFACE 3

#define VTX_TYPE_NAME "vertex_type"

unsigned getElmType(int numVtxPerElm) {
  if (numVtxPerElm == 4) {
    return apf::Mesh::TET;
  } else {
    fprintf(stderr, "unknown element type in %s\n", __func__);
    exit(EXIT_FAILURE);
  }
}

void readHeader(FILE* f, unsigned& nodes, unsigned& elms, unsigned& numVtxPerElm) {
  fscanf(f,"%d %d %d", &nodes, &elms, &numVtxPerElm);
}

void readCoords(FILE* f, unsigned numvtx, double* coordinates) {
  const double huge = 1024*1024;
  double min[3] = {huge,huge,huge};
  double max[3] = {-huge,-huge,-huge};
  for(unsigned i=0; i<numvtx; i++) {
    for(unsigned j=0; j<3; j++) {
      double pos = 0;
      fscanf(f, "%lf", &pos);
      coordinates[i*3+j] = pos;
      if( pos < min[j] )
        min[j] = pos;
      if( pos > max[j] )
        max[j] = pos;
    }
  }
  char d[3] = {'x','y','z'};
  fprintf(stderr, "boundary of the domain:\n");
  for(unsigned i=0; i<3; i++)
    fprintf(stderr, "%c %lf %lf \n", d[i], min[i], max[i]);
}

void readElements(FILE* f, unsigned numelms, int numVtxPerElm, int* elements) {
  unsigned i;
  std::map<int, int> count;
  for (i = 0; i < numelms*numVtxPerElm; i++) {
    int vtxid;
    fscanf(f, "%u", &vtxid);
    elements[i] = --vtxid; //export from matlab using 1-based indices
    count[elements[i]]++;
  }
  fprintf(stderr, "Number of vertices used %lu\n", count.size());
}

struct MeshInfo {
  double* coords;
  int* elements;
  unsigned elementType;
  unsigned numVerts;
  unsigned numElms;
  unsigned numVtxPerElm;
};

void readMesh(const char* meshfilename, MeshInfo& mesh) {
  FILE* f = fopen(meshfilename, "r"); 
  assert(f);
  readHeader(f,mesh.numVerts,mesh.numElms,mesh.numVtxPerElm);
  fprintf(stderr, "numVerts %u numElms %u numVtxPerElm %u\n", 
      mesh.numVerts, mesh.numElms, mesh.numVtxPerElm); 
  mesh.coords = new double[mesh.numVerts*3];
  readCoords(f, mesh.numVerts, mesh.coords);
  mesh.elements = new int [mesh.numElms*mesh.numVtxPerElm];
  readElements(f, mesh.numElms, mesh.numVtxPerElm, mesh.elements);
  mesh.elementType = getElmType(mesh.numVtxPerElm);
  fclose(f);
}

int* readVtx_type(const char* fname, unsigned numvtx) {
  FILE* f = fopen(fname, "r");
  assert(f);
  unsigned n;
  fscanf(f,"%u", &n);
  assert( n == numvtx );
  int* vtx_type_int = new int[numvtx];
  for(unsigned i = 0; i< numvtx; i++) {
    fscanf(f, "%d", &vtx_type_int[i]);
    assert(vtx_type_int[i] >= INTERIORTAG && vtx_type_int[i] <= TOP_PERIMETERTAG);
  }
  return vtx_type_int;
}

int main(int argc, char** argv)
{
  if( argc < 5 ) {
    printf("Usage: %s <GeomSim model .smd> <ascii mesh> <vertex classification> <output mesh> [<vertex field>...]\n", 
        argv[0]);
    return 0;
  }

  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  SimUtil_start();
  Sim_readLicenseFile(NULL);
  gmi_sim_start();
  gmi_register_mesh();
  gmi_register_sim();

  gmi_model* model = gmi_load(argv[1]);

  MeshInfo m;
  readMesh(argv[2],m);

  int* type = readVtx_type(argv[3], m.numVerts);
  delete [] type;

  const int dim = 3;
  apf::Mesh2* mesh = apf::makeEmptyMdsMesh(model, dim, false);
  apf::GlobalToVert outMap;
  apf::construct(mesh, m.elements, m.numElms, m.elementType, outMap);
  delete [] m.elements;
  apf::alignMdsRemotes(mesh);
  apf::deriveMdsModel(mesh);
  apf::setCoords(mesh, m.coords, m.numVerts, outMap);
  delete [] m.coords;
  outMap.clear();

  mesh->verify();
  mesh->writeNative(argv[3]);
  //apf::writeVtkFiles("after", mesh);

  mesh->destroyNative();
  apf::destroyMesh(mesh);
  gmi_destroy(model);
  gmi_sim_stop();
  Sim_unregisterAllKeys();
  SimUtil_stop();
  PCU_Comm_Free();
  MPI_Finalize();
}
