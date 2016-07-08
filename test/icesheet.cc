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

int* readVtxClassification(const char* fname, unsigned numvtx) {
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

apf::MeshTag* attachVtxClassification(apf::Mesh2* mesh, apf::GlobalToVert& vtxMap,
    unsigned numVerts, int* vtxClass) {
  apf::MeshTag* vtxTag = mesh->createIntTag(VTX_TYPE_NAME, 1);
  for(unsigned i=0; i<numVerts; i++)
    mesh->setIntTag(vtxMap[i], vtxTag, &(vtxClass[i]));
  delete [] vtxClass;
  return vtxTag;
}

bool isClassifiedOnBoundary(apf::Mesh2* mesh, apf::MeshEntity* face) {
  int numAdjElms = mesh->countUpward(face);
  assert( numAdjElms > 0 );
  if( numAdjElms == 2 )
    return false; // both regions exist -> interior
  else
    return true; // only one region -> exterior
}

/** \brief just get the first geometric model region - hack
*/
apf::ModelEntity* getMdlRgn(gmi_model* model) {
  gmi_iter* it = gmi_begin(model, 3);
  gmi_ent* rgn = gmi_next(model, it);
  assert(rgn);
  gmi_end(model, it);
  return (apf::ModelEntity*)rgn;
}

apf::ModelEntity* getMdlFace(apf::Mesh2* mesh, int tag) {
  apf::ModelEntity* face = mesh->findModelEntity(2,tag);
  assert(face);
  return face;
}

/** \brief brute force serach to find the intersection
 * \details models don't have that many edges bounding each face.... right???
 */
gmi_ent* getCommonEdge(gmi_set* a, gmi_set* b) {
  gmi_ent* match;
  unsigned found = 0;
  for(int i=0; i<a->n; i++) {
    for(int j=0; j<b->n; j++) {
      if( a->e[i] == b->e[j] ) {
        assert(!found);
        found++;
        match = a->e[i];
      }
    }
  }
  return match;
}

/** \brief get the model edge that bounds the model faces with the given tags
*/
apf::ModelEntity* getMdlEdge(apf::Mesh2* mesh, int faceTagA, int faceTagB) {
  gmi_model* model = mesh->getModel();
  gmi_ent* faceA = (gmi_ent*) getMdlFace(mesh,faceTagA);
  gmi_set* faceAedges = gmi_adjacent(model, faceA, 1);
  assert(faceAedges->n);
  gmi_ent* faceB = (gmi_ent*) getMdlFace(mesh,faceTagB);
  gmi_set* faceBedges = gmi_adjacent(model, faceB, 1);
  assert(faceBedges->n);
  gmi_ent* edge = getCommonEdge(faceAedges,faceBedges);
  assert(edge);
  gmi_free_set(faceAedges);
  gmi_free_set(faceBedges);
  return (apf::ModelEntity*) edge;
}

/* returns -1 for an interior face
 *          0 for a boundary face with at least one vertex with an interior tag
 *          1 for a marked boundary face
 */
int setModelClassification(gmi_model* model,
    apf::Mesh2* mesh,
    std::map<int,int> vtx_type_num,
    apf::MeshEntity* face,
    std::map<int,int>& faceClass) {
  if( ! isClassifiedOnBoundary(mesh,face) ) {
    mesh->setModelEntity(face,getMdlRgn(model));
    return -1;
  }
  //this is a boundary face
  assert(vtx_type_num.at(INTERIORTAG) == 0);
  //we only need to mark vertex that in on the boundary
  if(vtx_type_num[PERIMETERTAG]!=0) {
    //if perimeter point exist it's on the perimeter
    faceClass[PERIMETERFACE]++;
    mesh->setModelEntity(face,getMdlFace(mesh,PERIMETERFACE));
  }
  else {
    if((vtx_type_num[BOTTOMTAG] + vtx_type_num[BOTTOM_PERIMETERTAG]) == 3) {
      /*include vtx_type_num[BOTTOMTAG]=3;
        vtx_type_num[BOTTOMTAG]=2, vtx_type_num[BOTTOM_PERIMETERTAG])=1;
        vtx_type_num[BOTTOMTAG]=1, vtx_type_num[BOTTOM_PERIMETERTAG])=2;
        vtx_type_num[BOTTOM_PERIMETERTAG])=3;*/
      faceClass[BOTTOMFACE]++;
      mesh->setModelEntity(face,getMdlFace(mesh,BOTTOMFACE));
    }
    else {
      if((vtx_type_num[TOPTAG] + vtx_type_num[TOP_PERIMETERTAG]) == 3) {
        /*include vtx_type_num[TOPTAG]=3;
          vtx_type_num[TOPTAG]=2, vtx_type_num[TOP_PERIMETERTAG])=1;
          vtx_type_num[TOPTAG]=1, vtx_type_num[TOP_PERIMETERTAG])=2;
          vtx_type_num[BOTTOM_PERIMETERTAG])=3;*/
        faceClass[TOPFACE]++;
        mesh->setModelEntity(face,getMdlFace(mesh,TOPFACE));
      }
      else {
        faceClass[PERIMETERFACE]++;
        mesh->setModelEntity(face,getMdlFace(mesh,PERIMETERFACE));
      }
    }
  }
  return 1;
}

void setFaceClassification(gmi_model* model, apf::Mesh2* mesh, apf::MeshTag* vtxType) {
  int numbdryfaces = 0;
  int markedfaces = 0;
  int skippedfaces = 0;
  std::map<int,int> faceClass;

  apf::MeshIterator* it = mesh->begin(2);
  apf::MeshEntity* face;
  while( (face = mesh->iterate(it)) ) {
    apf::Downward verts;
    int n = mesh->getDownward(face, 0, verts);
    assert(n == 3);

    std::map<int, int> vtx_type_num;
    for(int i=INTERIORTAG; i<=TOP_PERIMETERTAG; i++)
      vtx_type_num[i] = 0;
    for(int i=0; i<3; i++){
      int value;
      mesh->getIntTag(verts[i], vtxType, &value);
      assert(value >= INTERIORTAG && value <= TOP_PERIMETERTAG);
      vtx_type_num[value]++;
    }
    int counttaggedvtx = 0;
    for(int i=INTERIORTAG; i<=TOP_PERIMETERTAG; i++)
      counttaggedvtx += vtx_type_num[i];
    assert(counttaggedvtx==3);
    int isSet = setModelClassification(model, mesh, vtx_type_num, face, faceClass);
    if( isSet == 1 )
      markedfaces++;
    if(isSet == 0)
      skippedfaces++;
    if( isSet == 0 || isSet == 1 )
      numbdryfaces++;
  }
  mesh->end(it);

  fprintf(stderr, "num boundary faces %d\n", numbdryfaces);
  fprintf(stderr, "num marked faces %d\n", markedfaces);
  fprintf(stderr, "num skipped faces %d\n", skippedfaces);
  int totmarkedfaces=0;
  for(int i=1; i<4; i++) {
    fprintf(stderr,"%d: %d\n", i, faceClass[i]);
    totmarkedfaces += faceClass[i];
  }
  assert(numbdryfaces == totmarkedfaces);
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

void getMeshEdgeTags(apf::Mesh2* mesh, apf::MeshEntity* edge,
    apf::MeshTag* t, int* tags) {
  apf::Downward verts;
  int n = mesh->getDownward(edge, 0, verts);
  assert(n == 2);

  for(int i=0; i<2; i++){
    int value;
    mesh->getIntTag(verts[i], t, &value);
    assert(value >= INTERIORTAG && value <= TOP_PERIMETERTAG);
    tags[i] = value;
  }
}

/** \brief set the mesh region classification
  \details hacked to set the classification to the same geometric model region
*/
void setEdgeClassification(gmi_model* model, apf::Mesh2* mesh, apf::MeshTag* t) {
  apf::ModelEntity* mdlRgn = getMdlRgn(model);
  apf::ModelEntity* mdlBotEdge = getMdlEdge(mesh,BOTTOMFACE,PERIMETERFACE);
  apf::ModelEntity* mdlTopEdge = getMdlEdge(mesh,TOPFACE,PERIMETERFACE);
  apf::ModelEntity* mdlTopFace = getMdlFace(mesh,TOPFACE);
  apf::ModelEntity* mdlBotFace = getMdlFace(mesh,BOTTOMFACE);
  apf::ModelEntity* mdlPerFace = getMdlFace(mesh,PERIMETERFACE);
  apf::MeshIterator* it = mesh->begin(1);
  apf::MeshEntity* edge;
  while( (edge = mesh->iterate(it)) ) {
    int tags[2];
    getMeshEdgeTags(mesh,edge,t,tags);
    // classified on model region
    if( tags[0] == INTERIORTAG || tags[1] == INTERIORTAG )
      mesh->setModelEntity(edge,mdlRgn);
    // classified on bottom perimeter
    else if( tags[0] == BOTTOM_PERIMETERTAG && tags[1] == BOTTOM_PERIMETERTAG )
      mesh->setModelEntity(edge,mdlBotEdge);
    // classified on top perimeter
    else if( tags[0] == TOP_PERIMETERTAG && tags[1] == TOP_PERIMETERTAG )
      mesh->setModelEntity(edge,mdlTopEdge);
    // classified on top face
    else if( tags[0] == TOPTAG || tags[1] == TOPTAG )
      mesh->setModelEntity(edge,mdlTopFace);
    // classified on bottom face
    else if( tags[0] == BOTTOMTAG || tags[1] == BOTTOMTAG )
      mesh->setModelEntity(edge,mdlBotFace);
    // classified on perimeter face
    else if( tags[0] == PERIMETERTAG || tags[1] == PERIMETERTAG )
      mesh->setModelEntity(edge,mdlPerFace);
    else {
      fprintf(stderr, "classification of edges fell through the conditional...exiting\n");
      exit(EXIT_FAILURE);
    }
  }
  mesh->end(it);
}

void setClassification(gmi_model* model, apf::Mesh2* mesh, apf::MeshTag* t) {
  setRgnClassification(model,mesh);
  setFaceClassification(model,mesh,t);
  setEdgeClassification(model,mesh,t);
  //setVtxClassification(model,mesh,t);
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


  const int dim = 3;
  apf::Mesh2* mesh = apf::makeEmptyMdsMesh(model, dim, false);
  apf::GlobalToVert outMap;
  apf::construct(mesh, m.elements, m.numElms, m.elementType, outMap);
  delete [] m.elements;
  apf::alignMdsRemotes(mesh);
  apf::deriveMdsModel(mesh);
  apf::setCoords(mesh, m.coords, m.numVerts, outMap);
  delete [] m.coords;
  int* vc = readVtxClassification(argv[3], m.numVerts);
  apf::MeshTag* vtxClass = attachVtxClassification(mesh, outMap, m.numVerts, vc);
  setClassification(model,mesh,vtxClass);
  outMap.clear();
  fprintf(stderr,"calling verify... flamesuit on\n");
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
