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

#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <regex>

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
//#define INTERIOR_REGION 0
int INTERIOR_REGION=0; // initialized but will be checked from read input

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
       INTERIOR_REGION=c-3000000;
       mesh->setModelEntity(v,getMdlRgn(model));
       //cint++;
    } else if (c >= 2000000) {
       mesh->setModelEntity(v,getMdlFace(mesh,c-2000000));
       //cface++;
    } else if (c >= 1000000) {
       mesh->setModelEntity(v,getMdlEdge(mesh,c-1000000));
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
  int k,ff;
  double distFromDebug1;
//  apf::Vector3 xd1(-0.54864, 7.44015e-06, 0.0397148 );
  apf::Vector3 xd1(-2.63644, -0.953687, 0.00762);
  apf::Vector3 xd2(0.914478, 0.0145401, 0.04 );
  apf::Vector3 dx1;
  apf::Vector3 dx2;
  apf::Vector3 tmp;
  apf::Vector3 Centroid;
  while( (e = mesh->iterate(it)) ) {
    //std::cout<<"Edge number "<<count++<<" with nodes"<<std::endl;
    mesh->getAdjacent(e, 0, verts);

    int nverts = verts.size();
    if(1==0) {
    Centroid=apf::getLinearCentroid(mesh,e);
    dx1=xd1-Centroid;
 //   dx2=xd2-Centroid;
    distFromDebug1=dx1[0]*dx1[0]
                  +dx1[1]*dx1[1]
                  +dx1[2]*dx1[2];
/*    distFromDebug2=dx2[0]*dx2[0]
                  +dx2[1]*dx2[1]
                  +dx2[2]*dx2[2];
*/
    }
    int cmax=-100;
    int cmin=100000000;
    for(int i=0; i<(int)verts.size(); i++) {
      mesh->getIntTag(verts[i],vtxClass,&c);
      cmax=std::max(cmax,c);
      cmin=std::min(cmin,c);
    }
    if(1==0 && distFromDebug1 < 1e-3) {
         fprintf(stderr, "%d %d %.15e %.15E %.15E \n", cmin, cmax, Centroid[0], Centroid[1], Centroid[2]);
         for (int i=0; i < nverts; i++) {
            mesh->getPoint(verts[i],0,tmp);
//            fprintf(stderr, "%d %.15e %.15E %.15E \n", i , tmp[0], tmp[1], tmp[2]);
         }
    }
    int dimMin=cmin/1000000;
    int dimMax=cmax/1000000;
    int tagMin=cmin-dimMin*1000000;
    int tagMax=cmax-dimMax*1000000;
    if (cmax >= 3000000) {
       INTERIOR_REGION=cmax-3000000;
       mesh->setModelEntity(e,getMdlRgn(model));
       //cint++;
    } else if (cmax >= 2000000) {  //max is a face
       if(cmin==cmax) { //min is same face -> correct face to cls
          mesh->setModelEntity(e,getMdlFace(mesh,cmax-2000000));
          //cface++;
       } else if (cmin >= 2000000) { // min is a DIFFERENT face -> interior
          mesh->setModelEntity(e,getMdlRgn(model));
       } else { 
//FAILS  ROLL OUR OWN         int res = gmi_is_in_closure_of(model,gmi_find(model,dimMin,tagMin), gmi_find(model,dimMax,tagMax));
          ff=-1;
          gmi_ent* ge=gmi_find(model,dimMin,tagMin);
          gmi_ent* gf =gmi_find(model,2,tagMax); // get the model face that goes with max
          gmi_set* Edges = gmi_adjacent(model,gf,1);
          k=0;
          while(k<((Edges->n)) && ff==-1){ // check all edges until one found
            if(dimMin==1) {
              if(ge==Edges->e[k]) ff=k;  // edges must be checked. 
            } else { // Verts probably can't fail but we still check
              gmi_set* Verts = gmi_adjacent(model,Edges->e[k],0);
              for (int j = 0; j < Verts->n; j++)  
                if(ge==Verts->e[j]) ff=j;
            }
            k++;
          }  
          if( ff!=-1 ) { // is it in cls
            mesh->setModelEntity(e,getMdlFace(mesh,cmax-2000000));
          } else  // edge not in closure so interior
          mesh->setModelEntity(e,getMdlRgn(model));
       }
    } else if (cmax >= 1000000) { // max is an edge
       if (cmin == cmax)  // cls on same edge 
          mesh->setModelEntity(e,getMdlEdge(mesh,cmax-1000000));
       else if (cmin >= 1000000) { // min is a different edge  and there is a face they must be in the closure of but which face is it 
         gmi_set* maxFaces = gmi_adjacent(model,gmi_find(model,1,tagMax),2);  
         gmi_set* minFaces = gmi_adjacent(model,gmi_find(model,1,tagMin),2);  
         for (int i = 0; i < maxFaces->n; i++) {
           for (int j = 0; j < minFaces->n; j++) {
             if(minFaces->e[i]==maxFaces->e[j]){
               int fftag=gmi_tag(model,maxFaces->e[j]);
               mesh->setModelEntity(e,getMdlFace(mesh,fftag));
             }
           }
         }
       } else mesh->setModelEntity(e,getMdlEdge(mesh,cmax-1000000)); // min is vtx thus max is correct edge to classify
    
    } else if (cmax < 1000000) { // two model verts so this is a 1 elm in z mesh
         gmi_set* maxEdges = gmi_adjacent(model,gmi_find(model,1,tagMax),1);  
         gmi_set* minEdges = gmi_adjacent(model,gmi_find(model,1,tagMin),1);  
         for (int i = 0; i < maxEdges->n; i++) {
           for (int j = 0; j < minEdges->n; j++) {
             if(minEdges->e[i]==maxEdges->e[j]){
               int fftag=gmi_tag(model,maxEdges->e[j]);
               mesh->setModelEntity(e,getMdlEdge(mesh,fftag));
             }
           }
         }
      
    } else { // should never get here since cmax < 10000000 is a vtx
            fprintf(stderr, "edge classification of these vert failed %d  %d \n", cmin, cmax);
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

  double distFromDebug1;
//  apf::Vector3 xd1(-0.54864, 7.44015e-06, 0.0397148 );
  apf::Vector3 xd1(-2.63644, -0.953687, 0.00762);
  apf::Vector3 xd2(0.914478, 0.0145401, 0.04 );
  apf::Vector3 dx1;
  apf::Vector3 dx2;
  apf::Vector3 tmp;
  apf::Vector3 Centroid;

  while( (f = mesh->iterate(it)) ) {
    mesh->getAdjacent(f, 0, verts);
    int nverts = verts.size();
    if(1==0) {
    Centroid=apf::getLinearCentroid(mesh,f);
    dx1=xd1-Centroid;
 //   dx2=xd2-Centroid;
    distFromDebug1=dx1[0]*dx1[0]
                  +dx1[1]*dx1[1]
                  +dx1[2]*dx1[2];
/*    distFromDebug2=dx2[0]*dx2[0]
                  +dx2[1]*dx2[1]
                  +dx2[2]*dx2[2];
*/
    }
    int cmin=100000000;
    int cmax=-100;
    int ctri[4];  // up to 4 points on a face
    for(int i=0; i<nverts; i++) {
      mesh->getIntTag(verts[i],vtxClass,&c);
      cmin=std::min(cmin,c);
      cmax=std::max(cmax,c);
      ctri[i]=c;
    }
//    if(std::min(distFromDebug1,distFromDebug2) < 1e-12) {
    if(1==0 && distFromDebug1 < 1e-3) {
         fprintf(stderr, "%d %d %.15e %.15E %.15E \n", cmin, cmax, Centroid[0], Centroid[1], Centroid[2]);
         for (int i=0; i < nverts; i++) {
            mesh->getPoint(verts[i],0,tmp);
//            fprintf(stderr, "%d %.15e %.15E %.15E \n", i , tmp[0], tmp[1], tmp[2]);
         }
    }
    if (cmax >= 3000000) { // at least one vertex is interior -> cls interior
       INTERIOR_REGION=cmax-3000000;
       mesh->setModelEntity(f,getMdlRgn(model));
       //cint++;
    } else if(cmin >= 2000000) { // all nodes on  model face(s?)
        if(cmax != cmin) { // all on faces but not all on same so classified on interior
          mesh->setModelEntity(f,getMdlRgn(model));
        } else { // all on same face so classify on that one
          mesh->setModelEntity(f,getMdlFace(mesh,cmax-2000000));
        }
    } else { // faces can ONLY be classified on model faces or interior but their vertices can be classified  on model faces, edge, or vertices (regions caught in if).  Consequently, the simplest logic is to loop over faces  and check if any face has all of this mesh face's verts model classification in its closure
      gmi_iter* gi=gmi_begin(model,2); // iterator over ALL the models faces
      gmi_ent* gf;
      gmi_ent* gt;
      int i,dimi,ff,tagi,k;
      int faceFound=0;
      int ifaceS=0;
      while ( (gf=gmi_next(model,gi)) && faceFound != nverts ) {
        faceFound=0;
        i=0;
        while(i<nverts && faceFound==i){ // check all verts but one not found == not in face
          ff=-1;
          dimi=ctri[i]/1000000;
          tagi=ctri[i]-dimi*1000000;
          gt=gmi_find(model,dimi,tagi); // model entity that vertex i is classified on
          // is_in_closure only works for simmetrix models so  roll our own
          if(dimi ==2 && gf==gt) ff=ifaceS;
          else if(dimi < 2) { // check this face's edges and those edges vertices
            gmi_ent* ge=gmi_find(model,dimi,tagi);
            gmi_set* Edges = gmi_adjacent(model,gf,1);
            k=0;
            while(k<((Edges->n)) && ff==-1){ // check all edges until one found
              if(dimi==1) {
                if(ge==Edges->e[k]) ff=k;  // edges must be checked. 
              } else { // Verts probably can't fail but we still check
                gmi_set* Verts = gmi_adjacent(model,Edges->e[k],0);
                for (int j = 0; j < Verts->n; j++)  
                   if(ge==Verts->e[j]) ff=j;
              }
              k++;
            }  
          }  
          if( ff!=-1 ) faceFound++;
          i++;
        }
        if(faceFound==nverts ) {
          int fftag=gmi_tag(model,gf);
          mesh->setModelEntity(f,getMdlFace(mesh,fftag));
        }
        ifaceS++;
      }
      gmi_end(model,gi);
      if(faceFound != nverts ) // none of the model face's closure held all verts classificaton  so interior
        mesh->setModelEntity(f,getMdlRgn(model));
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
  int maxmdlId=0;
  for(int i=0; i<localNumVtx; i++) {
    gmi_fscanf(f, 1, "%d",  &mdlId);
    (*classification)[i] = mdlId;
    maxmdlId=std::max(maxmdlId,mdlId); 
  }
  printf("maxmdlId,INTERIOR_REGION=%d %d\n",maxmdlId,INTERIOR_REGION);
  INTERIOR_REGION=maxmdlId-3000000; // correct if any vtx in region but will be negative if maxmdlId is on a face which is the case for 1 element in z... fix after reading model
  printf("maxmdlId,INTERIOR_REGION=%d %d\n",maxmdlId,INTERIOR_REGION);
}

void readFathers(FILE* f, int localNumVtx, int** fathers) {
  *fathers = new int[localNumVtx];
  rewind(f);
  int mdlId;
  for(int i=0; i<localNumVtx; i++) {
    gmi_fscanf(f, 1, "%d",  &mdlId);
    (*fathers)[i] = mdlId;
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

//static int starts_with(char const* with, char const* s) {
//  int lw;
//  int ls;
//  lw = strlen(with);
//  ls = strlen(s);
//  if (ls < lw)
//    return 0;
//  return strncmp(with, s, lw) == 0;
//}

bool seekPart(std::ifstream& f, const std::string& marker) {
  std::stringstream ss;
  ss << "^\\s+" << marker << "$";
  std::regex partId(ss.str());
  std::string line;
  while (std::getline(f, line)) {
    if (std::regex_match(line,partId)) {
      return true;
    }
  }
  return false;
}

struct BlockInfo {
  long numElms;
  int vtxPerElm;
};

std::vector<BlockInfo> readTopoBlockInfo(std::ifstream& f) {
  std::vector<BlockInfo> blocks;
  long blockSize;
  int vtxPerElement;

  std::string line;
  while (std::getline(f, line)) {
    std::istringstream iss(line);
    if (!(iss >> blockSize >> vtxPerElement)) { break; } // error
    blocks.push_back({blockSize,vtxPerElement});
  }
  return blocks;
}

void rewindStream(std::ifstream& f) {
  f.clear();
  f.seekg(0);
}

/**
fh = header file (there is only one for all processes), containing:
   Part0
       <numel_topo_1>   <NodesInElementTopo1>
       <nmel_topo_2>   < NodesInElementTopo2>
       ... for as many topos as are in  Part 0
   Repeat the above bock for each part.
**/
std::vector<BlockInfo> readHeader(std::ifstream& fh) {
  rewindStream(fh);
  const int self = PCU_Comm_Self();;
  bool ret = seekPart(fh, std::to_string(self));
  assert(ret);
  auto blockInfo = readTopoBlockInfo(fh);
  assert(blockInfo.size()>0);
  for(auto b : blockInfo) {
    std::cout << self << " " << b.numElms << " " << b.vtxPerElm << "\n";
  }
  return blockInfo;
}

/**
- f = part file, each part gets a file that contains a rectangular array, one for each topology
  present on that part, that provides element to vertexGlobalId connectivity in
  the order listed in the section of the header file for that part
**/
void readElements(std::ifstream& f, apf::Gid numElms,
    unsigned numVtxPerElm, apf::Gid* elements, bool rewind) {
  if(rewind) rewindStream(f);
  unsigned elmIdx = 0;
  apf::Gid* elmVtx = new apf::Gid[numVtxPerElm];
  for (int i = 0; i < numElms; i++) {
    for (unsigned j = 0; j < numVtxPerElm; j++)
      f >> elmVtx[j];
    for (unsigned j = 0; j < numVtxPerElm; j++) {
      const unsigned elmVtxIdx = elmIdx*numVtxPerElm+j;
      elements[elmVtxIdx] = --(elmVtx[j]); //export from matlab using 1-based indices
    }
    elmIdx++;
  }
  PCU_ALWAYS_ASSERT(numElms==elmIdx);
  delete [] elmVtx;
}

struct MeshInfo {
  double* coords;
  double* solution;
  std::vector<apf::Gid*> elements;
  apf::Gid* matches;
  int* classification;
  int* fathers2D;
  unsigned dim;
  std::vector<unsigned> elementType;
  apf::Gid numVerts;
  int localNumVerts;
  std::vector<apf::Gid> numElms;
  std::vector<unsigned> numVtxPerElm;
};

void readMesh(const char* meshfilename,
    const char* coordfilename,
    const char* matchfilename,
    const char* classfilename,
    const char* fathers2Dfilename,
    const char* solutionfilename,
    const char* connHeadfilename,
    MeshInfo& mesh) {

  mesh.dim = 3; //FIXME

  int self = PCU_Comm_Self();

  char filename[1024];
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
    readFathers(fff, mesh.localNumVerts, &(mesh.fathers2D)); // note we re-use classification reader
    fclose(fff);
  }

  if( strcmp(matchfilename, "NULL") ) {
    sprintf(filename, "%s.%d",matchfilename,self);
    FILE* fm = fopen(filename, "r");
    PCU_ALWAYS_ASSERT(fm);
    readMatches(fm, mesh.numVerts, mesh.localNumVerts, &(mesh.matches));
    fclose(fm);
  }

  std::stringstream ss;
  ss << meshfilename << "." << self;
  std::ifstream meshConnStream(ss.str());
  PCU_ALWAYS_ASSERT(meshConnStream.is_open());
  std::ifstream connHeadStream(connHeadfilename, std::ios::in);
  PCU_ALWAYS_ASSERT(connHeadStream.is_open());
  auto blockInfo = readHeader(connHeadStream);
  connHeadStream.close();
  bool rewind = true;
  for(auto b : blockInfo) {
    mesh.numElms.push_back(b.numElms);
    mesh.numVtxPerElm.push_back(b.vtxPerElm);
    apf::Gid* elements = new apf::Gid[b.numElms*b.vtxPerElm];
    readElements(meshConnStream, b.numElms, b.vtxPerElm, elements,rewind);
    rewind=false;
    mesh.elementType.push_back(getElmType(mesh.dim, b.vtxPerElm));
    mesh.elements.push_back(elements);
  }
  meshConnStream.close();
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
 //  1 element in z has no vertices to tell us the region ID so we have to find it from the tag number of the region that is adjacent to a face whose tag number will be the max when there is no region.  We know we have this problem when INTERIOR_REGION is negative
  if(INTERIOR_REGION < 0) {
    int maxmdlId=0;
    printf("maxmdlId,INTERIOR_REGION=%d %d\n",maxmdlId,INTERIOR_REGION);
    maxmdlId=1000000+INTERIOR_REGION;  // get back to maxmdlId found 
    printf("maxmdlId,INTERIOR_REGION=%d %d\n",maxmdlId,INTERIOR_REGION);
    gmi_set* Rgns=gmi_adjacent(model,gmi_find(model,2,maxmdlId),3);
    for(int i=0; i<((Rgns->n)); i++) 
      INTERIOR_REGION=gmi_tag(model,Rgns->e[i]);
  }
  printf("INTERIOR_REGION=%d\n",INTERIOR_REGION);
  apf::Mesh2* mesh = apf::makeEmptyMdsMesh(model, m.dim, isMatched);
  apf::GlobalToVert outMap;
  PCU_Debug_Open();
  for( size_t i=0; i< m.elements.size(); i++) {
    apf::assemble(mesh, m.elements[i], m.numElms[i], m.elementType[i], outMap);
    delete [] m.elements[i];
  }
  apf::finalise(mesh, outMap);
  apf::alignMdsRemotes(mesh);
  apf::deriveMdsModel(mesh);
  apf::setCoords(mesh, m.coords, m.localNumVerts, outMap);
  delete [] m.coords;
  if( isMatched ) {
    apf::setMatches(mesh, m.matches, m.localNumVerts, outMap);
    mesh->acceptChanges();
    delete [] m.matches;
  }
  apf::MeshTag* tc = setMappedTag(mesh, "classification", m.classification, 1,
      m.localNumVerts, outMap);
  { //debug - some verts don't have the classification tag...
    apf::MeshEntity* v;
    apf::MeshIterator* verts = mesh->begin(0);
    int i=0;
    while ((v = mesh->iterate(verts))) {
      if(!mesh->hasTag(v,tc)) {
        PCU_Debug_Print("%d missing tag\n", i);
      }
      i++;
    }
  }
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
  apf::writeVtkFiles("rendered",mesh);
  mesh->writeNative(argv[9]);
  if(noVerify != 1) mesh->verify();

  mesh->destroyNative();
  apf::destroyMesh(mesh);
  PCU_Comm_Free();
  MPI_Finalize();
}
