#include <gmi_mesh.h>
#include <gmi_null.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <apfConvert.h>
#include <apf.h>
#include <PCU.h>
#include <lionPrint.h>
#include <pcu_util.h>
#include <cstdlib>
#include <string.h>

unsigned getElmType(int numVtxPerElm) {
  if (numVtxPerElm == 4) {
    return apf::Mesh::TET;
  } else if (numVtxPerElm == 6) {
    return apf::Mesh::PRISM;
  } else if (numVtxPerElm == 8) {
    return apf::Mesh::HEX;
  } else {
    fprintf(stderr, "unknown element type in %s\n", __func__);
    exit(EXIT_FAILURE);
  }
}

bool skipLine(char* line) {
  // lines that start with either a '#' or a single white space
  // are skipped
  return (line[0] == '#' || line[0] == ' ' );
}

void getNumElms(FILE* f, unsigned& elms) {
  rewind(f);
  elms = 0;
  size_t linelimit = 1024;
  char* line = new char[linelimit];
  while( gmi_getline(&line,&linelimit,f) != -1 ) {
    if( ! skipLine(line) )
      elms++;
  }
  delete [] line;
}

void getNumVerts(FILE* f, unsigned& verts) {
  rewind(f);
  verts = 0;
  size_t linelimit = 1024;
  char* line = new char[linelimit];
  while( gmi_getline(&line,&linelimit,f) != -1 ) {
    if( ! skipLine(line) )
      verts++;
  }
  delete [] line;
}

void readCoords(FILE* f, unsigned numvtx, double* coordinates) {
  rewind(f);
  const double huge = 1024*1024;
  double min[3] = {huge,huge,huge};
  double max[3] = {-huge,-huge,-huge};
  for(unsigned i=0; i<numvtx; i++) {
    int id;
    gmi_fscanf(f, 1, "%d", &id);
    for(unsigned j=0; j<3; j++) {
      double pos = 0;
      gmi_fscanf(f, 1, "%lf", &pos);
      coordinates[i*3+j] = pos;
      if( pos < min[j] )
        min[j] = pos;
      if( pos > max[j] )
        max[j] = pos;
    }
  }
  char d[3] = {'x','y','z'};
  fprintf(stderr, "mesh bounding box:\n");
  for(unsigned i=0; i<3; i++)
    fprintf(stderr, "%c %lf %lf \n", d[i], min[i], max[i]);
}

void readMatches(FILE* f, unsigned numvtx, int* matches) {
  rewind(f);
  for(unsigned i=0; i<numvtx; i++) {
    int ignored, matchedVtx;
    gmi_fscanf(f, 2, "%d %d", &ignored, &matchedVtx); //export from matlab using 1-based indices
    PCU_ALWAYS_ASSERT( matchedVtx == -1 ||
        ( matchedVtx > 1 && matchedVtx <= static_cast<int>(numvtx) ));
    matches[i] = matchedVtx--;
  }
}

void readElements(FILE* f, unsigned numelms, unsigned numVtxPerElm,
    unsigned numVerts, int* elements) {
  rewind(f);
  unsigned i, j;
  std::map<int, int> count;
  for (i = 0; i < numelms; i++) {
    int elmid;
    gmi_fscanf(f, 1, "%u", &elmid);
    for (j = 0; j < numVtxPerElm; j++) {
      int elmVtxIdx = i*numVtxPerElm+j;
      int vtxid;
      gmi_fscanf(f, 1, "%u", &vtxid);
      elements[elmVtxIdx] = --vtxid; //export from matlab using 1-based indices
      count[elements[elmVtxIdx]]++;
    }
  }
  PCU_ALWAYS_ASSERT(count.size() == numVerts);
}

struct MeshInfo {
  double* coords;
  int* elements;
  int* matches;
  unsigned elementType;
  unsigned numVerts;
  unsigned numElms;
  unsigned numVtxPerElm;
};

void readMesh(const char* meshfilename,
    const char* coordfilename,
    const char* matchfilename,
    MeshInfo& mesh) {
  FILE* f = fopen(meshfilename, "r");
  PCU_ALWAYS_ASSERT(f);
  FILE* fc = fopen(coordfilename, "r");
  PCU_ALWAYS_ASSERT(fc);
  getNumElms(f,mesh.numElms);
  getNumVerts(fc,mesh.numVerts);
  fprintf(stderr, "numElms %u numVerts %u\n",
      mesh.numElms, mesh.numVerts);
  mesh.coords = new double[mesh.numVerts*3];
  readCoords(fc, mesh.numVerts, mesh.coords);
  fclose(fc);
  if( strcmp(matchfilename, "NULL") ) {
    FILE* fm = fopen(matchfilename, "r");
    PCU_ALWAYS_ASSERT(fm);
    mesh.matches = new int[mesh.numVerts];
    readMatches(fm, mesh.numVerts, mesh.matches);
    fclose(fm);
  }
  mesh.numVtxPerElm = 8; //hack!
  mesh.elements = new int [mesh.numElms*mesh.numVtxPerElm];
  readElements(f, mesh.numElms, mesh.numVtxPerElm, mesh.numVerts, mesh.elements);
  mesh.elementType = getElmType(mesh.numVtxPerElm);
  fclose(f);
}

void setMatches(apf::Mesh2* m, unsigned numVerts, int* matches,
    apf::GlobalToVert& globToVtx) {
  int self = PCU_Comm_Self();
  for(unsigned i=0; i<numVerts; i++) {
    if( matches[i] != -1 ) {
      m->addMatch(globToVtx[i], self, globToVtx[matches[i]]);
    }
  }
}

int main(int argc, char** argv)
{
  if( argc != 6 ) {
    printf("Usage: %s <ascii mesh connectivity .cnn> "
        "<ascii vertex coordinates .crd> "
        "<ascii vertex matching flag .match> "
        "<output model .dmg> <output mesh .smb>\n",
        argv[0]);
    return 0;
  }

  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  lion_set_verbosity(1);
  gmi_register_mesh();
  gmi_register_null();

  gmi_model* model = gmi_load(".null");

  MeshInfo m;
  readMesh(argv[1],argv[2],argv[3],m);

  bool isMatched = true;
  if( !strcmp(argv[3], "NULL") )
    isMatched = false;

  if(!PCU_Comm_Self())
    fprintf(stderr, "isMatched %d\n", isMatched);

  const int dim = 3;
  apf::Mesh2* mesh = apf::makeEmptyMdsMesh(model, dim, isMatched);
  apf::GlobalToVert outMap;
  apf::construct(mesh, m.elements, m.numElms, m.elementType, outMap);
  delete [] m.elements;
  apf::alignMdsRemotes(mesh);
  apf::deriveMdsModel(mesh);
  apf::setCoords(mesh, m.coords, m.numVerts, outMap);
  delete [] m.coords;
  if( isMatched ) {
    setMatches(mesh,m.numVerts,m.matches,outMap);
    delete [] m.matches;
  }
  outMap.clear();
  mesh->verify();

  gmi_write_dmg(model, argv[4]);
  mesh->writeNative(argv[5]);
  apf::writeVtkFiles("rendered",mesh);

  mesh->destroyNative();
  apf::destroyMesh(mesh);
  PCU_Comm_Free();
  MPI_Finalize();
}
