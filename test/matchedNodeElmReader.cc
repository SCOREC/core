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
#include <cassert>
#include <algorithm>
#include <apfBox.h>

/* from https://github.com/SCOREC/core/issues/205
0=fully interior of the volume
1-6 =classified on face (not edge or vertex)
11-22 = classified on model edge (not end points which are model vertices)
31-38 = classified on a model vertex.
*/

/* tags on vertices */
#define INTERIORTAG  0
#define FACE 1
#define FACE_LAST 6
#define EDGE 11
#define EDGE_LAST 22
#define VERTEX 31
#define VERTEX_LAST 38

/* model entity ids */
#define INTERIOR_REGION 0

apf::ModelEntity* getMdlRgn(gmi_model* model) {
  apf::ModelEntity* rgn = reinterpret_cast<apf::ModelEntity*>(
      gmi_find(model, 3, INTERIOR_REGION));
  PCU_ALWAYS_ASSERT(rgn);
  return rgn;
}

apf::ModelEntity* getMdlFace(apf::Mesh2* mesh, int tag) {
  apf::ModelEntity* face = mesh->findModelEntity(2,tag);
  PCU_ALWAYS_ASSERT(face);
  return face;
}

apf::ModelEntity* getMdlEdge(apf::Mesh2* mesh, int tag) {
  apf::ModelEntity* edge = mesh->findModelEntity(1,tag);
  PCU_ALWAYS_ASSERT(edge);
  return edge;
}

apf::ModelEntity* getMdlVtx(apf::Mesh2* mesh, int tag) {
  apf::ModelEntity* vertex = mesh->findModelEntity(0,tag);
  PCU_ALWAYS_ASSERT(vertex);
  return vertex;
}

void setVtxClassification(gmi_model* model, apf::Mesh2* mesh, apf::MeshTag* vtxClass) {
  (void)model;
  (void)mesh;
  (void)vtxClass;
  apf::MeshIterator* it = mesh->begin(0);
  apf::MeshEntity* v;
  //apf::Vector3 vCoord;
  int c;
  //int count=0,cint=0,cface=0,cedge=0,cvtx=0;
  while( (v = mesh->iterate(it)) ) {
    //mesh->getPoint(v, 0, vCoord);
    //std::cout<<"Coordinates: "<<vCoord[0]<<" "<<vCoord[1]<<" "<<vCoord[2]<<std::endl;
    mesh->getIntTag(v,vtxClass,&c);
    //std::cout<<"Returned tag is c= "<<c<<std::endl;
    //std::cout<<" "<<std::endl;
    //count++;
    if (c == INTERIORTAG) {
       mesh->setModelEntity(v,getMdlRgn(model));
       //cint++;
    } else if (c >= FACE && c <= FACE_LAST) {
       if (c == 1) { //face tag 1 corresponds to model face 0
          mesh->setModelEntity(v,getMdlFace(mesh,0));
       } else if (c == 2) { //face tag 2 corresponds to model face 1
          mesh->setModelEntity(v,getMdlFace(mesh,1));
       } else if (c == 3) { //face tag 3 corresponds to model face 3
          mesh->setModelEntity(v,getMdlFace(mesh,3));
       } else if (c == 4) { //face tag 4 corresponds to model face 4
          mesh->setModelEntity(v,getMdlFace(mesh,4));
       } else if (c == 5) { //face tag 5 corresponds to model face 2
          mesh->setModelEntity(v,getMdlFace(mesh,2));
       } else if (c == 6) { //face tag 6 corresponds to model face 5
          mesh->setModelEntity(v,getMdlFace(mesh,5));
       }
       //cface++;
    } else if (c >= EDGE && c <= EDGE_LAST) {
       if (c == 11) { //edge tag 11 corresponds to model edge 0
          mesh->setModelEntity(v,getMdlEdge(mesh,0));
       } else if (c == 12) {
          mesh->setModelEntity(v,getMdlEdge(mesh,2));
       } else if (c == 13) {
          mesh->setModelEntity(v,getMdlEdge(mesh,3));
       } else if (c == 14) {
          mesh->setModelEntity(v,getMdlEdge(mesh,1));
       } else if (c == 15) {
          mesh->setModelEntity(v,getMdlEdge(mesh,4));
       } else if (c == 16) {
          mesh->setModelEntity(v,getMdlEdge(mesh,5));
       } else if (c == 17) {
          mesh->setModelEntity(v,getMdlEdge(mesh,7));
       } else if (c == 18) {
          mesh->setModelEntity(v,getMdlEdge(mesh,6));
       } else if (c == 19) {
          mesh->setModelEntity(v,getMdlEdge(mesh,8));
       } else if (c == 20) {
          mesh->setModelEntity(v,getMdlEdge(mesh,10));
       } else if (c == 21) {
          mesh->setModelEntity(v,getMdlEdge(mesh,11));
       } else if (c == 22) {
          mesh->setModelEntity(v,getMdlEdge(mesh,9));
       }
       //cedge++;
    } else if (c >= VERTEX && c <= VERTEX_LAST) {
       if (c == 31) { //vertex tag 31 corresponds to model vertex 0
          mesh->setModelEntity(v,getMdlVtx(mesh,0));
       } else if (c == 32) { //vertex tag 32 corresponds to model vertex 1
          mesh->setModelEntity(v,getMdlVtx(mesh,1));
       } else if (c == 33) { //vertex tag 33 corresponds to model vertex 3
          mesh->setModelEntity(v,getMdlVtx(mesh,3));
       } else if (c == 34) { //vertex tag 34 corresponds to model vertex 2
          mesh->setModelEntity(v,getMdlVtx(mesh,2));
       } else if (c == 35) { //vertex tag 35 corresponds to model vertex 4
          mesh->setModelEntity(v,getMdlVtx(mesh,4));
       } else if (c == 36) { //vertex tag 36 corresponds to model vertex 5
          mesh->setModelEntity(v,getMdlVtx(mesh,5));
       } else if (c == 37) { //vertex tag 37 corresponds to model vertex 7
          mesh->setModelEntity(v,getMdlVtx(mesh,7));
       } else if (c == 38) { //vertex tag 38 corresponds to model vertex 6
          mesh->setModelEntity(v,getMdlVtx(mesh,6));
       }
       //cvtx++;
    }
  }
  //std::cout<<"count is "<<cvtx<<std::endl;
  mesh->end(it);
}

void setEdgeClassification(gmi_model* model, apf::Mesh2* mesh,apf::MeshTag* vtxClass) {
  (void)model;
  (void)mesh;
  (void)vtxClass;
  apf::MeshIterator* it = mesh->begin(1);
  apf::MeshEntity* e;
  //apf::Vector3 vCoord;
  int c;
  //int count=0;
  apf::Adjacent verts;
  while( (e = mesh->iterate(it)) ) {
    //std::cout<<"Edge number "<<count++<<" with nodes"<<std::endl;
    mesh->getAdjacent(e, 0, verts);
    int cmin=100;
    for(size_t i=0; i<verts.size(); i++) {
      mesh->getIntTag(verts[i],vtxClass,&c);
      //mesh->getPoint(verts[i], 0, vCoord);
      //std::cout<<vCoord[0]<<" "<<vCoord[1]<<" "<<vCoord[2]<<std::endl;
      cmin=std::min(cmin,c);
    }
    //std::cout<<"has classification "<<cmin<<std::endl;
    //std::cout<<" "<<std::endl;
    if (cmin == INTERIORTAG) {
       mesh->setModelEntity(e,getMdlRgn(model));
       //cint++;
    } else if (cmin >= FACE && cmin <= FACE_LAST) {
       if (cmin == 1) { //face tag 1 corresponds to model face 0
          mesh->setModelEntity(e,getMdlFace(mesh,0));
       } else if (cmin == 2) { //face tag 2 corresponds to model face 1
          mesh->setModelEntity(e,getMdlFace(mesh,1));
       } else if (cmin == 3) { //face tag 3 corresponds to model face 3
          mesh->setModelEntity(e,getMdlFace(mesh,3));
       } else if (cmin == 4) { //face tag 4 corresponds to model face 4
          mesh->setModelEntity(e,getMdlFace(mesh,4));
       } else if (cmin == 5) { //face tag 5 corresponds to model face 2
          mesh->setModelEntity(e,getMdlFace(mesh,2));
       } else if (cmin == 6) { //face tag 6 corresponds to model face 5
          mesh->setModelEntity(e,getMdlFace(mesh,5));
       }
       //cface++;
    } else if (cmin >= EDGE && cmin <= EDGE_LAST) {
       if (cmin == 11) { //edge tag 11 corresponds to model edge 0
          mesh->setModelEntity(e,getMdlEdge(mesh,0));
       } else if (cmin == 12) {
          mesh->setModelEntity(e,getMdlEdge(mesh,2));
       } else if (cmin == 13) {
          mesh->setModelEntity(e,getMdlEdge(mesh,3));
       } else if (cmin == 14) {
          mesh->setModelEntity(e,getMdlEdge(mesh,1));
       } else if (cmin == 15) {
          mesh->setModelEntity(e,getMdlEdge(mesh,4));
       } else if (cmin == 16) {
          mesh->setModelEntity(e,getMdlEdge(mesh,5));
       } else if (cmin == 17) {
          mesh->setModelEntity(e,getMdlEdge(mesh,7));
       } else if (cmin == 18) {
          mesh->setModelEntity(e,getMdlEdge(mesh,6));
       } else if (cmin == 19) {
          mesh->setModelEntity(e,getMdlEdge(mesh,8));
       } else if (cmin == 20) {
          mesh->setModelEntity(e,getMdlEdge(mesh,10));
       } else if (cmin == 21) {
          mesh->setModelEntity(e,getMdlEdge(mesh,11));
       } else if (cmin == 22) {
          mesh->setModelEntity(e,getMdlEdge(mesh,9));
       }
       //cedge++;
    } else if (cmin >= VERTEX && cmin <= VERTEX_LAST) {
       if (cmin == 31) { //vertex tag 31 corresponds to model vertex 0
          mesh->setModelEntity(e,getMdlVtx(mesh,0));
       } else if (cmin == 32) { //vertex tag 32 corresponds to model vertex 1
          mesh->setModelEntity(e,getMdlVtx(mesh,1));
       } else if (cmin == 33) { //vertex tag 33 corresponds to model vertex 3
          mesh->setModelEntity(e,getMdlVtx(mesh,3));
       } else if (cmin == 34) { //vertex tag 34 corresponds to model vertex 2
          mesh->setModelEntity(e,getMdlVtx(mesh,2));
       } else if (cmin == 35) { //vertex tag 35 corresponds to model vertex 4
          mesh->setModelEntity(e,getMdlVtx(mesh,4));
       } else if (cmin == 36) { //vertex tag 36 corresponds to model vertex 5
          mesh->setModelEntity(e,getMdlVtx(mesh,5));
       } else if (cmin == 37) { //vertex tag 37 corresponds to model vertex 7
          mesh->setModelEntity(e,getMdlVtx(mesh,7));
       } else if (cmin == 38) { //vertex tag 38 corresponds to model vertex 6
          mesh->setModelEntity(e,getMdlVtx(mesh,6));
       }
       //cvtx++;
    }
  }
  mesh->end(it);
}

/* if any of four vertices are classified on region -> region
 * else on model face and it is impossible to have more than one face in the 4
 * vertices classification
 * */
void setFaceClassification(gmi_model* model, apf::Mesh2* mesh, apf::MeshTag* vtxClass) {
  (void)model;
  (void)mesh;
  (void)vtxClass;

  apf::MeshIterator* it = mesh->begin(2);
  apf::MeshEntity* f;
  int c;
/* What was here before was commented when we went to simple min rule
  apf::Adjacent verts;
  while( (f = mesh->iterate(it)) ) {
    m->getAdjacent(f, 0, verts) = 0;
    bool hasRgClass = false;
    for(int i=0; i<verts.size(); i++) {
      m->getIntTag(i,vtxClass,&c);
      if( c == INTERIORTAG )
        mesh->setModelEntity(f,mdlRgn);
      }
      
    }
  }
*/
  apf::Adjacent verts;
  while( (f = mesh->iterate(it)) ) {
    mesh->getAdjacent(f, 0, verts);
    int cmin=100;
    for(size_t i=0; i<verts.size(); i++) {
      mesh->getIntTag(verts[i],vtxClass,&c);
      cmin=std::min(cmin,c);
    }
    if (cmin == INTERIORTAG) {
       mesh->setModelEntity(f,getMdlRgn(model));
       //cint++;
    } else if (cmin >= FACE && cmin <= FACE_LAST) {
       if (cmin == 1) { //face tag 1 corresponds to model face 0
          mesh->setModelEntity(f,getMdlFace(mesh,0));
       } else if (cmin == 2) { //face tag 2 corresponds to model face 1
          mesh->setModelEntity(f,getMdlFace(mesh,1));
       } else if (cmin == 3) { //face tag 3 corresponds to model face 3
          mesh->setModelEntity(f,getMdlFace(mesh,3));
       } else if (cmin == 4) { //face tag 4 corresponds to model face 4
          mesh->setModelEntity(f,getMdlFace(mesh,4));
       } else if (cmin == 5) { //face tag 5 corresponds to model face 2
          mesh->setModelEntity(f,getMdlFace(mesh,2));
       } else if (cmin == 6) { //face tag 6 corresponds to model face 5
          mesh->setModelEntity(f,getMdlFace(mesh,5));
       }
       //cface++;
    } else if (cmin >= EDGE && cmin <= EDGE_LAST) {
       if (cmin == 11) { //edge tag 11 corresponds to model edge 0
          mesh->setModelEntity(f,getMdlEdge(mesh,0));
       } else if (cmin == 12) {
          mesh->setModelEntity(f,getMdlEdge(mesh,2));
       } else if (cmin == 13) {
          mesh->setModelEntity(f,getMdlEdge(mesh,3));
       } else if (cmin == 14) {
          mesh->setModelEntity(f,getMdlEdge(mesh,1));
       } else if (cmin == 15) {
          mesh->setModelEntity(f,getMdlEdge(mesh,4));
       } else if (cmin == 16) {
          mesh->setModelEntity(f,getMdlEdge(mesh,5));
       } else if (cmin == 17) {
          mesh->setModelEntity(f,getMdlEdge(mesh,7));
       } else if (cmin == 18) {
          mesh->setModelEntity(f,getMdlEdge(mesh,6));
       } else if (cmin == 19) {
          mesh->setModelEntity(f,getMdlEdge(mesh,8));
       } else if (cmin == 20) {
          mesh->setModelEntity(f,getMdlEdge(mesh,10));
       } else if (cmin == 21) {
          mesh->setModelEntity(f,getMdlEdge(mesh,11));
       } else if (cmin == 22) {
          mesh->setModelEntity(f,getMdlEdge(mesh,9));
       }
       //cedge++;
    } else if (cmin >= VERTEX && cmin <= VERTEX_LAST) {
       if (cmin == 31) { //vertex tag 31 corresponds to model vertex 0
          mesh->setModelEntity(f,getMdlVtx(mesh,0));
       } else if (cmin == 32) { //vertex tag 32 corresponds to model vertex 1
          mesh->setModelEntity(f,getMdlVtx(mesh,1));
       } else if (cmin == 33) { //vertex tag 33 corresponds to model vertex 3
          mesh->setModelEntity(f,getMdlVtx(mesh,3));
       } else if (cmin == 34) { //vertex tag 34 corresponds to model vertex 2
          mesh->setModelEntity(f,getMdlVtx(mesh,2));
       } else if (cmin == 35) { //vertex tag 35 corresponds to model vertex 4
          mesh->setModelEntity(f,getMdlVtx(mesh,4));
       } else if (cmin == 36) { //vertex tag 36 corresponds to model vertex 5
          mesh->setModelEntity(f,getMdlVtx(mesh,5));
       } else if (cmin == 37) { //vertex tag 37 corresponds to model vertex 7
          mesh->setModelEntity(f,getMdlVtx(mesh,7));
       } else if (cmin == 38) { //vertex tag 38 corresponds to model vertex 6
          mesh->setModelEntity(f,getMdlVtx(mesh,6));
       }
       //cvtx++;
    }
  }
  mesh->end(it);
}

/** \brief set the mesh region classification
  \details hacked to set the classification to the same geometric model region
*/
void setRgnClassification(gmi_model* model, apf::Mesh2* mesh) {
  apf::ModelEntity* mdlRgn = getMdlRgn(model);
  apf::MeshIterator* it = mesh->begin(3);
  apf::MeshEntity* rgn;
  while( (rgn = mesh->iterate(it)) )
    mesh->setModelEntity(rgn,mdlRgn);
  mesh->end(it);
}

void setClassification(gmi_model* model, apf::Mesh2* mesh, apf::MeshTag* t) {
  setRgnClassification(model,mesh);
  setFaceClassification(model,mesh,t);
  setEdgeClassification(model,mesh,t);
  setVtxClassification(model,mesh,t);
  mesh->acceptChanges();
}

void getLocalRange(unsigned total, unsigned& local,
    long& first, long& last) {
  const int self = PCU_Comm_Self();
  const int peers = PCU_Comm_Peers();
  local = total/peers;
  if( self == peers-1 ) //last rank
    if( local*peers < total )
      local += total - local*peers;
  first = PCU_Exscan_Long(local);
  last = first+local;
}

void printElmTypeError(int dim, int numVtxPerElm) {
  fprintf(stderr, "unknown element type for"
      "dim %d and numVtxPerElm %d in %s\n",
      dim, numVtxPerElm, __func__);
}

unsigned getElmType(int dim, int numVtxPerElm) {
  if (dim == 2) {
    if (numVtxPerElm == 3)
      return apf::Mesh::TRIANGLE;
    if (numVtxPerElm == 4)
      return apf::Mesh::QUAD;
    else {
      printElmTypeError(dim, numVtxPerElm);
      exit(EXIT_FAILURE);
    }
  } else if (dim == 3) {
    if (numVtxPerElm == 4)
      return apf::Mesh::TET;
    else if (numVtxPerElm == 6)
      return apf::Mesh::PRISM;
    else if (numVtxPerElm == 8)
      return apf::Mesh::HEX;
    else {
      printElmTypeError(dim, numVtxPerElm);
      exit(EXIT_FAILURE);
    }
  } else {
    printElmTypeError(dim, numVtxPerElm);
    exit(EXIT_FAILURE);
  }
}

bool skipLine(char* line) {
  // lines that start with either a '#' or a single white space
  // are skipped
  return (line[0] == '#' || line[0] == ' ' );
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

void readClassification(FILE* f, unsigned numVtx, int** classification) {
  long firstVtx, lastVtx;
  unsigned localNumVtx;
  getLocalRange(numVtx,localNumVtx,firstVtx,lastVtx);
  *classification = new int[localNumVtx];
  rewind(f);
  int vidx = 0;
  for(unsigned i=0; i<numVtx; i++) {
    int id;
    int mdlId;
    gmi_fscanf(f, 2, "%d %d", &id, &mdlId);
    //std::cout<<"read id is "<<id<<std::endl;
    //std::cout<<"read model id is "<<mdlId<<std::endl;
    if( i >= firstVtx && i < lastVtx ) {
      (*classification)[vidx] = mdlId;
      vidx++;
    }
  }
  /*std::cout<<"numVtx is "<<numVtx<<std::endl;
  for (int i=0;i<numVtx;i++) {
      std::cout<<"vtx num "<<i<<" has class "<<(*classification)[i]<<std::endl;
  }*/
}

void readCoords(FILE* f, unsigned numvtx, unsigned& localnumvtx, double** coordinates) {
  long firstVtx, lastVtx;
  getLocalRange(numvtx, localnumvtx,firstVtx,lastVtx);
  *coordinates = new double[localnumvtx*3];
  rewind(f);
  int vidx = 0;
  for(unsigned i=0; i<numvtx; i++) {
    int id;
    double pos[3];
    gmi_fscanf(f, 4, "%d %lf %lf %lf", &id, pos+0, pos+1, pos+2);
    if( i >= firstVtx && i < lastVtx ) {
      for(unsigned j=0; j<3; j++)
        (*coordinates)[vidx*3+j] = pos[j];
      vidx++;
    }
  }
}

void readMatches(FILE* f, unsigned numvtx, int** matches) {
  long firstVtx, lastVtx;
  unsigned localnumvtx;
  getLocalRange(numvtx, localnumvtx, firstVtx, lastVtx);
  fprintf(stderr, "%d readMatches numvtx %d localnumvtx %u firstVtx %ld lastVtx %ld\n",
      PCU_Comm_Self(), numvtx, localnumvtx, firstVtx, lastVtx);
  *matches = new int[localnumvtx];
  rewind(f);
  int vidx = 0;
  int gid, matchedVtx;
  int i = 0;
  while( 2 == fscanf(f, "%d %d", &gid, &matchedVtx) ) {
    if( i >= firstVtx && i < lastVtx ) {
      PCU_ALWAYS_ASSERT( matchedVtx == -1 ||
          ( matchedVtx >= 1 && matchedVtx <= static_cast<int>(numvtx) ));
      if( matchedVtx != -1 )
        --matchedVtx;
      if( matchedVtx == 66350 || matchedVtx == 65075 ) {
        fprintf(stderr, "%d reader found match %d at gid %d i %d vidx %d\n",
            PCU_Comm_Self(), matchedVtx, gid, i, vidx);
      }
      (*matches)[vidx] = matchedVtx;
      vidx++;
    }
    i++;
  }
}

void readElements(FILE* f, unsigned &dim, unsigned& numElms,
    unsigned& numVtxPerElm, unsigned& localNumElms, int** elements) {
  rewind(f);
  int dimHeader[2];
  gmi_fscanf(f, 2, "%u %u", dimHeader, dimHeader+1);
  assert( dimHeader[0] == 1 && dimHeader[1] == 1);
  gmi_fscanf(f, 1, "%u", &dim);
  gmi_fscanf(f, 2, "%u %u", &numElms, &numVtxPerElm);
  long firstElm, lastElm;
  getLocalRange(numElms, localNumElms, firstElm, lastElm);
  *elements = new int[localNumElms*numVtxPerElm];
  unsigned i, j;
  unsigned elmIdx = 0;
  int* elmVtx = new int[numVtxPerElm];
  for (i = 0; i < numElms; i++) {
    int ignored;
    gmi_fscanf(f, 1, "%u", &ignored);
    for (j = 0; j < numVtxPerElm; j++)
      gmi_fscanf(f, 1, "%u", elmVtx+j);
    if (i >= firstElm && i < lastElm) {
      for (j = 0; j < numVtxPerElm; j++) {
        const unsigned elmVtxIdx = elmIdx*numVtxPerElm+j;
        (*elements)[elmVtxIdx] = --(elmVtx[j]); //export from matlab using 1-based indices
      }
      elmIdx++;
    }
  }
  delete [] elmVtx;
}

struct MeshInfo {
  double* coords;
  int* elements;
  int* matches;
  int* classification;
  int* fathers2D;
  unsigned dim;
  unsigned elementType;
  unsigned numVerts;
  unsigned localNumVerts;
  unsigned numElms;
  unsigned localNumElms;
  unsigned numVtxPerElm;
};

void readMesh(const char* meshfilename,
    const char* coordfilename,
    const char* matchfilename,
    const char* classfilename,
    const char* fathers2Dfilename,
    MeshInfo& mesh) {
  FILE* fc = fopen(coordfilename, "r");
  PCU_ALWAYS_ASSERT(fc);
  getNumVerts(fc,mesh.numVerts);
  if(!PCU_Comm_Self())
    fprintf(stderr, "numVerts %u\n", mesh.numVerts);
  readCoords(fc, mesh.numVerts, mesh.localNumVerts, &(mesh.coords));
  fclose(fc);

  FILE* ff = fopen(classfilename, "r");
  PCU_ALWAYS_ASSERT(ff);
  readClassification(ff, mesh.numVerts, &(mesh.classification));
  fclose(ff);

  //add an argument to readMesh for the fathers2D
  FILE* fff = fopen(fathers2Dfilename, "r");
  PCU_ALWAYS_ASSERT(fff);
  readClassification(fff, mesh.numVerts, &(mesh.fathers2D)); // note we re-use classification reader
  fclose(fff);

  if( strcmp(matchfilename, "NULL") ) {
    FILE* fm = fopen(matchfilename, "r");
    PCU_ALWAYS_ASSERT(fm);
    readMatches(fm, mesh.numVerts, &(mesh.matches));
    fclose(fm);
  }

  FILE* f = fopen(meshfilename, "r");
  PCU_ALWAYS_ASSERT(f);
  readElements(f, mesh.dim, mesh.numElms, mesh.numVtxPerElm,
      mesh.localNumElms, &(mesh.elements));
  mesh.elementType = getElmType(mesh.dim, mesh.numVtxPerElm);
  fclose(f);
}


int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  lion_set_verbosity(1);
  if( argc != 8 ) {
    if( !PCU_Comm_Self() ) {
      printf("Usage: %s <ascii mesh connectivity .cnn> "
          "<ascii vertex coordinates .crd> "
          "<ascii vertex matching flag .match> "
          "<ascii vertex classification flag .class> "
          "<ascii vertex fathers2D flag .fathers2D> "
          "<output model .dmg> <output mesh .smb>\n",
          argv[0]);
    }
    return 0;
  }

  gmi_register_mesh();
  gmi_register_null();


  double t0 = PCU_Time();
  MeshInfo m;
  readMesh(argv[1],argv[2],argv[3],argv[4],argv[5],m);

  bool isMatched = true;
  if( !strcmp(argv[3], "NULL") )
    isMatched = false;

  if(!PCU_Comm_Self())
    fprintf(stderr, "isMatched %d\n", isMatched);

  //gmi_model* model = gmi_load(".null");
  gmi_model* model = apf::makeMdsBox(2,2,2,1,1,1,0);
  apf::Mesh2* mesh = apf::makeEmptyMdsMesh(model, m.dim, isMatched);
  apf::GlobalToVert outMap;
  apf::construct(mesh, m.elements, m.localNumElms, m.elementType, outMap);
  delete [] m.elements;
  apf::alignMdsRemotes(mesh);
  apf::deriveMdsModel(mesh);
  /*for (int i=0; i<81; i++) {
  std::cout<<m.coords[i]<<std::endl;
  }*/
  apf::setCoords(mesh, m.coords, m.localNumVerts, outMap);
  delete [] m.coords;
  if( isMatched ) {
    apf::setMatches(mesh, m.matches, m.localNumVerts, outMap);
    mesh->acceptChanges();
    delete [] m.matches;
  }
  apf::MeshTag* tc = setIntTag(mesh, "classification", m.classification, 1,
      m.localNumVerts, outMap);
  setClassification(model,mesh,tc);
  apf::removeTagFromDimension(mesh, tc, 0);
  mesh->destroyTag(tc);
 
  apf::MeshTag* tf = setIntTag(mesh, "fathers2D", m.fathers2D, 1,
      m.localNumVerts, outMap);
  apf::removeTagFromDimension(mesh, tf, 0);
  mesh->destroyTag(tf);

  if(!PCU_Comm_Self())
    fprintf(stderr, "seconds to create mesh %.3f\n", PCU_Time()-t0);
  mesh->verify();

  outMap.clear();
  gmi_write_dmg(model, argv[6]);
  mesh->writeNative(argv[7]);
  apf::writeVtkFiles("rendered",mesh);

  mesh->destroyNative();
  apf::destroyMesh(mesh);
  PCU_Comm_Free();
  MPI_Finalize();
}
