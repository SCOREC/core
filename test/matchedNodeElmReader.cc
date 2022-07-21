#include <gmi_mesh.h>
#include <gmi_null.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <apfConvert.h>
#include <apfConvertTags.h>
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
//
// New approach:  mesher will add 0 to mdlConvert's dmg tag number for vertices,
//                          1000000 to mdlConvert's dmg tag number for edges,
//                          2000000 to mdlConvert's dmg tag number for faces,
//                          3000000 to mdlConvert's dmg tag number for regions
    if (c >= 3000000 ) {
       mesh->setModelEntity(v,getMdlRgn(model));
       //cint++;
    } else if (c >= 2000000) {
       mesh->setModelEntity(v,getMdlFace(mesh,c-2000000));
       //cface++;
    } else if (c >= 1000000) {
       mesh->setModelEntity(v,getMdlEdge(mesh,c-100000));
       //cedge++;
    } else {
       mesh->setModelEntity(v,getMdlVtx(mesh,c));
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
    int cmax=-100;
    int cmin=100000000;
    for(size_t i=0; i<verts.size(); i++) {
      mesh->getIntTag(verts[i],vtxClass,&c);
      //mesh->getPoint(verts[i], 0, vCoord);
      //std::cout<<vCoord[0]<<" "<<vCoord[1]<<" "<<vCoord[2]<<std::endl;
      cmax=std::max(cmax,c);
      cmin=std::min(cmin,c);
    }
    //std::cout<<"has classification "<<cmin<<std::endl;
    //std::cout<<" "<<std::endl;
    if (cmax >= 3000000) {
       mesh->setModelEntity(e,getMdlRgn(model));
       //cint++;
    } else if (cmax >= 2000000) {  //max is a face
       if(cmin==cmax) { //min is same face -> correct face to cls
          mesh->setModelEntity(e,getMdlFace(mesh,cmax-2000000));
          //cface++;
       } else if (cmin >= 2000000) { // min is a DIFFERENT face -> interior
          mesh->setModelEntity(e,getMdlRgn(model));
       } else if (cmin >= 1000000) { // max is face min is an edge 
         if( GF_inClosure(getMdlFace(mesh,cmax), getMdlEdge(mesh,cmin)) ) { // is it in cls
           mesh->setModelEntity(e,getMdlFace(mesh,cmax-2000000));
         } else  // edge not in closure so interior
           mesh->setModelEntity(e,getMdlRgn(model));
       } else if( GF_inClosure(getMdlFace(mesh,cmax), getMdlVtx(mesh,cmin)) ) { // second must be vtx but is that vtx in closure  ACTUALLY IT MUST BE SO UNNECESSARY
            mesh->setModelEntity(e,getMdlFace(mesh,cmax-2000000));
       } else 
           mesh->setModelEntity(e,getMdlRgn(model));
    } else if (cmax >= 1000000) { // max is an edge
       if (cmin == cmax)  // cls on same edge 
          mesh->setModelEntity(e,getMdlEdge(mesh,cmax-1000000));
       else if (cmin >= 1000000) { // min is a different edge  and there is a face they must be in the closure of but which face is it 
         pPList maxFaces = GE_faces(getMdlEdge(mesh,cmax-1000000);  
         pPList minFaces = GE_faces(getMdlEdge(mesh,cmin-1000000);  
         for (int i = 0; i < PList_size(maxFaces); i++) {
           pGFace maxFace = (pGFace) PList_item(maxFaces, i);
           for (int j = 0; j < PList_size(minFaces); j++) {
             pGFace minFace = (pGFace) PList_item(minFaces, j);
             if(minFace==maxFace) mesh->setModelEntity(e,minFace));
           }
         }
       } else mesh->setModelEntity(e,getMdlEdge(mesh,cmax-1000000)); // min is vtx thus max is correct edge to classify
    } else { // should never get her
       std:cout <<"classification of edge on a vertex is not valid";
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
  apf::Adjacent verts;
  while( (f = mesh->iterate(it)) ) { // loop over all mesh faces
    mesh->getAdjacent(f, 0, verts); 
    size_t nverts = verts.size();
    int cmin=100000000;
    int cmax=-100;
    int ctri[4];  // up to 4 points on a face
    for(size_t i=0; i<nverts; i++) {
      mesh->getIntTag(verts[i],vtxClass,&c);
      cmin=std::min(cmin,c);
      cmax=std::max(cmax,c);
      ctri[i]=c;
    }
    if (cmax >= 3000000) { // at least one vertex is interior -> cls interior
       mesh->setModelEntity(e,getMdlRgn(model));
       //cint++;
    } else if(cmin >= 2000000) { // all nodes on  model face(s?)
        if(cmax != cmin) { // all on faces but not all on same so classified on interior
          mesh->setModelEntity(f,getMdlRgn(model));
        } else { // all on same face so classify on that one
          mesh->setModelEntity(f,getMdlFace(mesh,cmax-2000000));
        }
    } else { // faces can ONLY be classified on model faces or interior but their vertices can be classified  on model faces, edge, or vertices (regions caught in if).  Consequently, the simplest logic is to loop over faces  and check  if all verts in closure
      int faceFound=0;
      pGFace gfaceFound;
      GFIter gfIter=GM_faceIter(model);
      while ( faceFound != verts.size() && (gface=GFIter_next(gfIter))) {
        faceFound=0;
        for(size_t i=0; i< nverts; i++) { // check if each vert is in the cls of faceh
          if(ctri[i] >= 2000000 {   // i is a face 
             if( GF_inClosure(gface, getMdlFace(mesh,ctri[i]))) faceFound++;
	  } else if(ctri[i] >= 1000000 ) { // i is an edge 
             if( GF_inClosure(gface, getMdlEdge(mesh,ctri[i]))) faceFound++;
          } else {// i is a vertex
             if( GF_inClosure(gface, getMdlVtx(mesh,ctri[i]))) faceFound++;
          }
        }
        if(faceFound==verts.size() ) {
          gfaceFound=gface;
          mesh->setModelEntity(f,gfaceFound);
          // does C++ have a concept of exiting the wile loop
        }
      }         
     }
     if(faceFound != verts.size() ) 
            fprintf(stderr, "face classification of these vert classification failed %d %d  %d \n", cmin, cmid, cmax);
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



void getLocalRange(apf::Gid total, int& local,
    apf::Gid& first, apf::Gid& last) {
  const int self = PCU_Comm_Self();
  const int peers = PCU_Comm_Peers();
  local = total/peers; 
  if( self == peers-1 ) {  //last rank
    apf::Gid lp=local*peers;
    if( lp < total ){
      apf::Gid lpd;
      lpd= total - lp;
      local += lpd;
    }
  }
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

void getNumVerts(FILE* f, apf::Gid& verts) {
  rewind(f);
  gmi_fscanf(f, 1, "%ld",  &verts);
}

void readClassification(FILE* f, int localNumVtx, int** classification) {
  *classification = new int[localNumVtx];
  rewind(f);
  int mdlId;
  for(int i=0; i<localNumVtx; i++) {
    gmi_fscanf(f, 1, "%d",  &mdlId);
    (*classification)[i] = mdlId;
  }
}

void readCoords(FILE* f, int& localnumvtx, double** coordinates) {
  *coordinates = new double[localnumvtx*3];
  for(int i=0; i<localnumvtx; i++) {
    double pos[3];
    gmi_fscanf(f, 3, "%lf %lf %lf", pos+0, pos+1, pos+2);
    for(unsigned j=0; j<3; j++)
        (*coordinates)[i*3+j] = pos[j];
  }
}

void readSolution(FILE* f, int& localnumvtx, double** solution) {
  *solution = new double[localnumvtx*5];
  rewind(f);
  for(int i=0; i<localnumvtx; i++) {
    double pos[5];
    pos[4]=0; //temperature
    gmi_fscanf(f, 4, "%lf %lf %lf %lf", pos+0, pos+1, pos+2, pos+3);
    for(unsigned j=0; j<5; j++)
      (*solution)[i*5+j] = pos[j];
  }
}

void readMatches(FILE* f, apf::Gid numvtx, int localnumvtx, apf::Gid** matches) {
  fprintf(stderr, "%d readMatches numvtx %ld localnumvtx %d \n",
      PCU_Comm_Self(), numvtx, localnumvtx);
  *matches = new apf::Gid[localnumvtx];
  rewind(f);
  apf::Gid matchedVtx;
  for(int i=0; i<localnumvtx; i++) {
    gmi_fscanf(f, 1, "%ld",  &matchedVtx);
    PCU_ALWAYS_ASSERT( matchedVtx == -1 ||
       ( matchedVtx >= 1 && matchedVtx <= numvtx ));
    if( matchedVtx != -1 )
        --matchedVtx;
    (*matches)[i] = matchedVtx;
  }
// I think the above will perform better than the code commented out below
//  int vidx = 0;
//  while( 1 == fscanf(f, "%ld", &matchedVtx) ) {
//    PCU_ALWAYS_ASSERT( matchedVtx == -1 ||
//       ( matchedVtx >= 1 && matchedVtx <= numvtx ));
//    if( matchedVtx != -1 )
//        --matchedVtx;
//    (*matches)[vidx] = matchedVtx;
//    vidx++;
//  }
}

void readElements(FILE* f, FILE* fh, unsigned &dim,  apf::Gid& numElms,
    unsigned& numVtxPerElm, int& localNumElms, apf::Gid** elements) {
  rewind(f);
  rewind(fh);
  int dimHeader[2];
  gmi_fscanf(fh, 2, "%u %u", dimHeader, dimHeader+1);
//  assert( dimHeader[0] == 1 && dimHeader[1] == 1);
  gmi_fscanf(fh, 1, "%u", &dim);
  gmi_fscanf(fh, 2, "%ld %u", &numElms, &numVtxPerElm);
  int self = PCU_Comm_Self();;
  for (int j=0; j< self+1;j++)
     gmi_fscanf(fh, 2, "%d %u", &localNumElms, &numVtxPerElm);
  *elements = new apf::Gid[localNumElms*numVtxPerElm];
  int i;
  unsigned j;
  unsigned elmIdx = 0;
  apf::Gid* elmVtx = new apf::Gid[numVtxPerElm];
  for (i = 0; i < localNumElms; i++) {
    for (j = 0; j < numVtxPerElm; j++)
      gmi_fscanf(f, 1, "%ld", elmVtx+j);
    for (j = 0; j < numVtxPerElm; j++) {
      const unsigned elmVtxIdx = elmIdx*numVtxPerElm+j;
      (*elements)[elmVtxIdx] = --(elmVtx[j]); //export from matlab using 1-based indices
    }
    elmIdx++;
  }
  delete [] elmVtx;
}

struct MeshInfo {
  double* coords;
  double* solution;
  apf::Gid* elements;
  apf::Gid* matches;
  int* classification;
  int* fathers2D;
  unsigned dim;
  unsigned elementType;
  apf::Gid numVerts;
  int localNumVerts;
  apf::Gid numElms;
  int localNumElms;
  unsigned numVtxPerElm;
};

void readMesh(const char* meshfilename,
    const char* coordfilename,
    const char* matchfilename,
    const char* classfilename,
    const char* fathers2Dfilename,
    const char* solutionfilename,
    const char* connHeadfilename,
    MeshInfo& mesh) {

  int self = PCU_Comm_Self();

  char filename[64];
  sprintf(filename, "%s.%d",coordfilename,self);
    
  FILE* fc = fopen(filename , "r");
  PCU_ALWAYS_ASSERT(fc);
  getNumVerts(fc,mesh.numVerts);
  mesh.localNumVerts=mesh.numVerts;
  mesh.numVerts=PCU_Add_Long(mesh.numVerts);
  
  if(!PCU_Comm_Self())
    fprintf(stderr, "numVerts %ld\n", mesh.numVerts);
  readCoords(fc, mesh.localNumVerts, &(mesh.coords));
  fclose(fc);
 
  if(0==1) {
  sprintf(filename, "%s.%d",solutionfilename,self);
  FILE* fs = fopen(filename, "r");
  PCU_ALWAYS_ASSERT(fs);
  readSolution(fs, mesh.localNumVerts, &(mesh.solution));
  fclose(fs);
  }

  sprintf(filename, "%s.%d",classfilename,self);
  FILE* ff = fopen(filename, "r");
  PCU_ALWAYS_ASSERT(ff);
  readClassification(ff, mesh.localNumVerts, &(mesh.classification));
  fclose(ff);

  if( strcmp(fathers2Dfilename, "NULL") ) {
  //add an argument to readMesh for the fathers2D
    sprintf(filename, "%s.%d",fathers2Dfilename,self);
    FILE* fff = fopen(filename, "r");
    PCU_ALWAYS_ASSERT(fff);
    readClassification(fff, mesh.localNumVerts, &(mesh.fathers2D)); // note we re-use classification reader
    fclose(fff);
  }

  if( strcmp(matchfilename, "NULL") ) {
    sprintf(filename, "%s.%d",matchfilename,self);
    FILE* fm = fopen(filename, "r");
    PCU_ALWAYS_ASSERT(fm);
    readMatches(fm, mesh.numVerts, mesh.localNumVerts, &(mesh.matches));
    fclose(fm);
  }

  sprintf(filename, "%s.%d",meshfilename,self);
  FILE* f = fopen(filename, "r");
  FILE* fh = fopen(connHeadfilename, "r");
  PCU_ALWAYS_ASSERT(f);
  PCU_ALWAYS_ASSERT(fh);
  readElements(f,fh, mesh.dim, mesh.numElms, mesh.numVtxPerElm,
      mesh.localNumElms, &(mesh.elements));
  mesh.elementType = getElmType(mesh.dim, mesh.numVtxPerElm);
  fclose(f);
}


int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  lion_set_verbosity(1);
  int noVerify=0;    // maintain default of verifying if not explicitly requesting it off
  if( argc < 10 ) {
    if( !PCU_Comm_Self() ) {
      printf("Usage: %s <ascii mesh connectivity .cnn> "
          "<ascii vertex coordinates .crd> "
          "<ascii vertex matching flag .match> "
          "<ascii vertex classification flag .class> "
          "<ascii vertex fathers2D flag .fathers2D> "
          "<ascii solution flag .soln> "
          "<ascii conn header> "
          "<output model .dmg> <output mesh .smb>"
          "turn off verify mesh if equal 1 (on if you give nothing)\n",
          argv[0]);
    }
    return 0;
  }

  gmi_register_mesh();
  gmi_register_null();

  if( argc == 11 ) noVerify=atoi(argv[10]);

  double t0 = PCU_Time();
  MeshInfo m;
  readMesh(argv[1],argv[2],argv[3],argv[4],argv[5],argv[6],argv[7],m);

  bool isMatched = true;
  if( !strcmp(argv[3], "NULL") )
    isMatched = false;

  if(!PCU_Comm_Self())
    fprintf(stderr, "isMatched %d\n", isMatched);

  //gmi_model* model = gmi_load(".null");
  gmi_model* model = gmi_load("outModelr.dmg");
//  gmi_model* model = apf::makeMdsBox(2,2,2,1,1,1,0);
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
  apf::MeshTag* tc = setMappedTag(mesh, "classification", m.classification, 1,
      m.localNumVerts, outMap);
  setClassification(model,mesh,tc);
  apf::removeTagFromDimension(mesh, tc, 0);
  mesh->destroyTag(tc);
 
  if( strcmp(argv[5], "NULL") ) {
    apf::MeshTag* tf = setMappedTag(mesh, "fathers2D", m.fathers2D, 1,
      m.localNumVerts, outMap);
    (void) tf;
  } else if(!PCU_Comm_Self())
    fprintf(stderr, "fathers2D not requested \n");

  //mesh->destroyTag(tf);

  if(0==1) {
  apf::MeshTag* ts = setMappedTag(mesh, "solution", m.solution, 5,
      m.localNumVerts, outMap);
  (void) ts;
  }


  /* // Print the father2D tags
  apf::MeshEntity* v;
  apf::MeshIterator* it = mesh->begin(0);
  apf::MeshTag* t = mesh->findTag("fathers2D");
  if (t==NULL) {std::cout<<"Didn't find tag"<<std::endl;}
  int tagNum; 
  int count = 0;
  while ((v = mesh->iterate(it))) { // loop over mesh vertices
    mesh->getIntTag(v,t,&tagNum);
    std::cout<<"Tag number "<<tagNum<<std::endl;
    count++;
  }
  mesh->end(it);*/

  if(!PCU_Comm_Self())
    fprintf(stderr, "seconds to create mesh %.3f\n", PCU_Time()-t0);
  if(noVerify != 1) mesh->verify();

  outMap.clear();
  gmi_write_dmg(model, argv[8]);
  mesh->writeNative(argv[9]);
  if(noVerify != 0) mesh->verify();
  apf::writeVtkFiles("rendered",mesh);

  mesh->destroyNative();
  apf::destroyMesh(mesh);
  PCU_Comm_Free();
  MPI_Finalize();
}
