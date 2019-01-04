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
  unsigned elementType;
  unsigned numVerts;
  unsigned numElms;
  unsigned numVtxPerElm;
};

void readMesh(const char* meshfilename, const char* coordfilename, MeshInfo& mesh) {
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
  mesh.numVtxPerElm = 8; //hack!
  mesh.elements = new int [mesh.numElms*mesh.numVtxPerElm];
  readElements(f, mesh.numElms, mesh.numVtxPerElm, mesh.numVerts, mesh.elements);
  mesh.elementType = getElmType(mesh.numVtxPerElm);
  fclose(f);
}

int* readIntArray(const char* fname, unsigned len) {
  FILE* f = fopen(fname, "r");
  PCU_ALWAYS_ASSERT(f);
  unsigned n;
  gmi_fscanf(f, 1, "%u", &n);
  PCU_ALWAYS_ASSERT( n == len );
  int* data = new int[len];
  for(unsigned i = 0; i< len; i++)
    gmi_fscanf(f, 1, "%d", &data[i]);
  fclose(f);
  return data;
}

double* readScalarArray(const char* fname, unsigned len) {
  FILE* f = fopen(fname, "r");
  PCU_ALWAYS_ASSERT(f);
  unsigned n;
  gmi_fscanf(f, 1, "%u", &n);
  PCU_ALWAYS_ASSERT( n == len );
  double* data = new double[len];
  for(unsigned i = 0; i< len; i++)
    gmi_fscanf(f, 1, "%lf", &data[i]);
  fclose(f);
  return data;
}

std::string getName(const char* fname) {
  using std::string;
  std::string in(fname);
  if( in.find("vertex_type") != string::npos ) {
    return string("vertex_type");
  } else if( in.find("basal_friction") != string::npos ) {
    return string("basal_friction");
  } else if( in.find("temperature") != string::npos ) {
    return string("temperature");
  } else if( in.find("surface_elevation") != string::npos ) {
    return string("surface_height");
  } else if( in.find("solution_x") != string::npos ) {
    return string("solution_x");
  } else if( in.find("solution_y") != string::npos ) {
    return string("solution_y");
  } else {
    fprintf(stderr, "unknown field name in %s\n", __func__);
    exit(EXIT_FAILURE);
  }
}

apf::MeshTag* attachVtxTag(apf::Mesh2* mesh, const char* fname,
    apf::GlobalToVert& vtxMap) {
  unsigned numverts = mesh->count(0);
  std::string fldname = getName(fname);
  fprintf(stderr, "attaching %s\n", fldname.c_str());
  int* arr = readIntArray(fname, numverts);
  apf::MeshTag* tag = mesh->createIntTag(fldname.c_str(), 1);
  for(unsigned i=0; i<numverts; i++)
    mesh->setIntTag(vtxMap[i], tag, &(arr[i]));
  delete [] arr;
  return tag;
}

/** \brief attach a field to the vertices
 *  \detail all the fields being attached need to exist with the mesh so the
 *  field pointer is not returned
 */
void attachVtxField(apf::Mesh2* mesh, const char* fname,
    apf::GlobalToVert& vtxMap) {
  std::string fldname = getName(fname);
  unsigned numverts = mesh->count(0);
  fprintf(stderr, "attaching %s\n", fldname.c_str());
  double* arr = readScalarArray(fname, numverts);
  apf::Field* fld = apf::createFieldOn(mesh,fldname.c_str(),apf::SCALAR);
  for(unsigned i=0; i<numverts; i++)
    apf::setScalar(fld,vtxMap[i],0,arr[i]);
  delete [] arr;
}

void mergeSolutionFields(apf::Mesh2* mesh) {
  apf::Field* x = mesh->findField("solution_x");
  apf::Field* y = mesh->findField("solution_y");
  apf::Field* xy = apf::createFieldOn(mesh,"Solution",apf::VECTOR);
  apf::MeshIterator* it = mesh->begin(0);
  apf::MeshEntity* vtx;
  while( (vtx = mesh->iterate(it)) ) {
    apf::Vector3 v;
    v[0] = apf::getScalar(x,vtx,0);
    v[1] = apf::getScalar(y,vtx,0);
    v[2] = 0.0;
    apf::setVector(xy, vtx, 0, v);
  }
  mesh->end(it);
  apf::destroyField(x);
  apf::destroyField(y);
}

void detachVtxTag(apf::Mesh2* mesh, apf::MeshTag* t) {
  apf::MeshIterator* it = mesh->begin(0);
  apf::MeshEntity* vtx;
  while( (vtx = mesh->iterate(it)) )
    mesh->removeTag(vtx,t);
  mesh->end(it);
  mesh->destroyTag(t);
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

  //gmi_model* model = gmi_load(".null");

  MeshInfo m;
  readMesh(argv[1],argv[2],m);

  /*
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

  gmi_write_dmg(model, argv[4]);
  mesh->writeNative(argv[5]);
  apf::writeVtkFiles("rendered",mesh);

  mesh->destroyNative();
  apf::destroyMesh(mesh);
  */
  PCU_Comm_Free();
  MPI_Finalize();
}
