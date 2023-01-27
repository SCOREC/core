#include <PCU.h>
#include <lionPrint.h>
#include <MeshSim.h>
#include <SimPartitionedMesh.h>
#include "SimParasolidKrnl.h"
#include <SimAdvMeshing.h>
#include <SimUtil.h>
#include <apfSIM.h>
#include <apfMDS.h>
#include <gmi.h>
#include <gmi_sim.h>
#include <apf.h>
#include <apfConvert.h>
#include <apfMesh2.h>
#include <apfNumbering.h>
#include <apfShape.h>
#include <ma.h>
#include <pcu_util.h>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <cassert>
#include <getopt.h>
#include <string.h>
#include <stdio.h>

using namespace std;




apf::Field* convert_my_tag(apf::Mesh* m, apf::MeshTag* t) {
  apf::MeshEntity* vtx;
  apf::MeshIterator* it = m->begin(0);
  apf::Field* f = apf::createFieldOn(m, "fathers2D_field", apf::SCALAR);
  int vals[1];
  double vals_d;
  while ((vtx = m->iterate(it))) {
    m->getIntTag(vtx, t, vals);
    vals_d = vals[0];
    apf::setScalar(f, vtx, 0, vals_d);
  }
  m->end(it);
  return f;
}

static void attachOrder(apf::Mesh* m)
{
  apf::numberOverlapDimension(m, "sim_order", m->getDimension());
}

static void fixMatches(apf::Mesh2* m)
{
  if (m->hasMatching()) {
    if (apf::alignMdsMatches(m))
      printf("fixed misaligned matches\n");
    else
      printf("matches were aligned\n");
    PCU_ALWAYS_ASSERT( ! apf::alignMdsMatches(m));
  }
}

static void fixPyramids(apf::Mesh2* m)
{
  if (m->getDimension() != 3)
    return; /* no pyramids exist in 2D */
  if (apf::countEntitiesOfType(m, apf::Mesh::HEX))
    return; /* meshadapt can't even look at hexes */
  ma::Input* in = ma::makeAdvanced(ma::configureIdentity(m));
  in->shouldCleanupLayer = true;
  ma::adapt(in);
}

const char* gmi_path = NULL;
const char* gmi_native_path = NULL;
const char* sms_path = NULL;
const char* smb_path = NULL;
int should_log = 0;
int should_fix_pyramids = 1;
int should_attach_order = 0;
const char* extruRootPath = NULL;
int ExtruRootId =0;
bool found_bad_arg = false;

void getConfig(int argc, char** argv) {

  opterr = 0;

  static struct option long_opts[] = {
    {"no-pyramid-fix", no_argument, &should_fix_pyramids, 0},
    {"attach-order", no_argument, &should_attach_order, 1},
    {"enable-log", no_argument, &should_log, 2},
    {"model-face-root", required_argument, 0, 'e'},
    {"native-model", required_argument, 0, 'n'},
    {0, 0, 0, 0}  // terminate the option array
  };

  const char* usage=""
    "[options] <model file> <simmetrix mesh> <scorec mesh>\n"
    "options:\n"
    "  --no-pyramid-fix                Disable quad-connected pyramid tetrahedronization\n"
    "  --attach-order                  Attach the Simmetrix element order as a Numbering\n"
    "  --enable-log                    Enable Simmetrix logging\n"
    "  --model-face-root=/path/to/file ASCII input file with one integer per line listing the face ids that are the roots of mesh extrusions from SimModeler\n"
    "  --native-model=/path/to/model   Load the native Parasolid or ACIS model that the GeomSim model uses\n";

  int option_index = 0;
  while(1) {
    int c = getopt_long(argc, argv, "", long_opts, &option_index);
    if (c == -1) break; //end of options
    switch (c) {
      case 0: // pyramid fix flag
      case 1: // attach order flag
      case 2: // enable simmetrix logging
        break;
      case 'e':
        extruRootPath = optarg;
        break;
      case 'n':
        gmi_native_path = optarg;
        break;
      case '?':
        if (!PCU_Comm_Self())
          printf ("warning: skipping unrecognized option\n");
        break;
      default:
        if (!PCU_Comm_Self())
          printf("Usage %s %s", argv[0], usage);
        exit(EXIT_FAILURE);
    }
  }

  if(argc-optind != 3) {
    if (!PCU_Comm_Self())
      printf("Usage %s %s", argv[0], usage);
    exit(EXIT_FAILURE);
  }
  int i=optind;
  gmi_path = argv[i++];
  sms_path = argv[i++];
  smb_path = argv[i++];
//  gmi_native_path= "/nobackup/uncompressed/Models/GustWing/2dCTS/SlicedGeomtryV12/sms2mds/geom.xmt_txt";
  gmi_native_path= "geom.xmt_txt";
  if (!PCU_Comm_Self()) {
    printf ("fix_pyramids %d attach_order %d enable_log %d extruRootPath %s\n",
            should_fix_pyramids, should_attach_order, should_log, extruRootPath);
    printf ("native-model \'%s\' model \'%s\' simmetrix mesh \'%s\' output mesh \'%s\'\n",
      gmi_native_path, gmi_path, sms_path, smb_path);
  }
}

// put the extrude tagging here which 1) loops over the mesh faces classified on the model face that is the root of the extrude
// create a tag on vertices fathers
// get the list of mesh rootfaces classified on the source geometric model face
// for each srcFace in rootfaces
// get the ids of downward adjacent vertices, store that as an array of size 3
// get the upward adjacent region srcRgn
// call Extrusion_3DRegionsAndLayerFaces(srcRgn,...)
// for each face in the returned list of faces
// get the downward adjacent vertices of face - they will be in the same order as the srcFace ids
// set the fathers tag
// assert that the x,y coordinates of each vertex matches the srcFace vertex coordinates within some relaxed
// tolerance - sanity check my assumption that face-to-vtx adjaceny is always the same order
void addFathersTag(pGModel simModel, pParMesh sim_mesh, apf::Mesh* simApfMesh, const char* extrusionFaceFile) {
  if(!extrusionFaceFile) return;
  // create a tag on vertices fathers
  pMeshDataId myFather = MD_newMeshDataId( "fathers2D");
  
  pPList listV,listVn,faces,regions;
  pFace face;
  pRegion region;
  pVertex vrts[4];
  int dir, err;
  int count2D=0;
  pGRegion gregion;
  pGFace gface;
  pGEdge gedge;
  pGFace ExtruRootFace=NULL;
  pVertex entV;
  pMesh meshP= PM_mesh	(sim_mesh, 0 );	

  char coordfilename[64];
  char cnnfilename[64];
  sprintf(coordfilename, "geom.crd");
  sprintf(cnnfilename, "geom.cnn");
  FILE* fcr = fopen(coordfilename, "w");
  FILE* fcn = fopen(cnnfilename, "w");

  FILE* fid = fopen(extrusionFaceFile, "r"); // helper file that contains all faces with extrusions
  assert(fid);
  double VdisTol=1e-12;
  double EdisTol=1e-7;
  double FdisTol=1e-7;
  while(1 == fscanf(fid,"%d",&ExtruRootId)) {
    fprintf(stderr,"ExtruRootId= %d \n",ExtruRootId);
    //find the root face of the extrusion
    GFIter gfIter=GM_faceIter(simModel);
    while ( (gface=GFIter_next(gfIter))) {
      int id = GEN_tag(gface);
      if(id==ExtruRootId) ExtruRootFace=gface;
    }
    assert(ExtruRootFace != NULL);

    FIter fIter = M_classifiedFaceIter( meshP, ExtruRootFace, 0 ); // 0 says I don't want closure
    while ((face = FIter_next(fIter))) {
      dir=1;
      listV= F_vertices(face, dir);
      void *iter = 0;        // Must initialize to 0
      int i=0;
      while ((entV =(pVertex)PList_next(listV, &iter))) { //loop over plist of vertices
        // Process each item in list
        vrts[i] = (pVertex)entV;
        i++;
      }
      int nvert=i;
      PList_delete(listV);

      double coordNewPt[nvert][3];
      for(i=0; i< nvert ; i++) {
        int* markedData;
        if(!EN_getDataPtr((pEntity)vrts[i],myFather,(void**)&markedData)){  // not sure about marked yet
          gType  vClassDim;
          pGEntity vConG;
          int foundESTag = 0;
          int foundETag = 0;
          int foundEETag = 0;
          pPList gRegions,gFaces,gEdges;
          double coordGVSelf[3];
          double coordGVOther[3];
          double coordGVOther1[3];
          double dx,dy;
          vClassDim=V_whatInType(vrts[i]);
          vConG=V_whatIn(vrts[i]);
          int echeck;
          double parFace[2]; 
          double aLow[2];
          double aHigh[2];
          double pLow, pHigh;
          switch(vClassDim) {
            case 0:   // classified on vert so vert->edge->vert
              foundESTag = GEN_tag( (pGVertex) vConG ); // found Extrusion Start Tag
              gEdges = GV_edges( (pGVertex) vConG ); // pPList of model edges adjacent to root model vertex
              for(int j = 0; j < PList_size( gEdges ); j++ ){
                pGEdge gEdge = (pGEdge) PList_item( gEdges , j ); // candidate edge
                pGVertex gVert0 = GE_vertex( gEdge , 0 );
                pGVertex gVert1 = GE_vertex( gEdge , 1 );
                if( gVert0 == (pGVertex) vConG) { //1 is at other end of edge
                   GV_point( gVert1 , coordGVOther );
                   GV_point( gVert0 , coordGVSelf );
                   dx = coordGVOther[0] - coordGVSelf[0];
                   dy = coordGVOther[1] - coordGVSelf[1];
                   if( dx*dx + dy*dy < VdisTol ) {
                      foundETag = GEN_tag( gEdge );
                      foundEETag = GEN_tag( gVert1 );
                   }
                } else { // 0 is at the other edge
                   GV_point( gVert0 , coordGVOther );
                   GV_point( gVert1 , coordGVSelf );
                   dx = coordGVOther[0] - coordGVSelf[0];
                   dy = coordGVOther[1] - coordGVSelf[1];
                   if( dx*dx + dy*dy < VdisTol) {
                      foundETag = GEN_tag(gEdge);
                      foundEETag = GEN_tag(gVert0);
                   }
                }
              } 
              PList_delete(gEdges);
              break; 
            case 1:   // classified on edge so edge->face->edge
              foundESTag = GEN_tag( (pGEdge) vConG ); // found Extrusion Start Tag
              GE_parRange ( (pGEdge) vConG, &pLow, &pHigh);
              parFace[0] = ( pLow + pHigh ) * 0.5;
              echeck = GE_point( (pGEdge) vConG , parFace[0], coordGVSelf ); 
              gFaces = GE_faces( (pGEdge) vConG); // pPList of model faces adjacent to root model edge
              for(int j = 0; j < PList_size( gFaces ); j++ ){
                pGFace gFace = (pGFace) PList_item( gFaces , j ); // candidate face
                gEdges = GF_edges( gFace ); // pPList of model edges of jth adjacent face
                for(int k = 0; k < PList_size( gEdges ); k++ ){ // loop over that pPlist
                  pGEdge gEdge = (pGEdge)  PList_item( gEdges , k ); // candidate edge on candidate face
                  if( gEdge != (pGEdge) vConG ) { // exclude root classified edge
                    GE_parRange ( gEdge, &pLow, &pHigh);
                    parFace[0] = ( pLow + pHigh ) * 0.5;
                    echeck = GE_point( gEdge , parFace[0], coordGVOther ); 
                    dx = coordGVOther[0] - coordGVSelf[0];
                    dy = coordGVOther[1] - coordGVSelf[1];
                    if( dx*dx + dy*dy < EdisTol) {
                       foundETag = GEN_tag( gFace ); // found Extruded Tag
                       foundEETag = GEN_tag( gEdge );   // found Extrusion End Tag
                    }
                  }
                }
                PList_delete(gEdges);
              }
              PList_delete(gFaces);

              break; 
            case 2:   // classified on face so face->region->face
              foundESTag = GEN_tag( (pGFace) vConG ); // found Extrusion Start Tag
              GF_parRange ( (pGFace) vConG, 0, &pLow, &pHigh);
              parFace[0] = ( pLow + pHigh ) * 0.5;
              GF_parRange ( (pGFace) vConG, 1, &pLow, &pHigh);
              parFace[1] = ( pLow + pHigh ) * 0.5;
              echeck = GF_point( (pGFace) vConG, parFace , coordGVSelf );
              gRegions = GF_regions( (pGFace) vConG ); //pPList of model regions adjacent to root model Fac
              pGRegion gRegion = (pGRegion) PList_item( gRegions , 0 ); // there can be only one for extrusions
              int  currentTagRegion = GEN_tag( gRegion ); // candidate region's tag (debugging only)
              gFaces = GR_faces( gRegion );
              for( int j = 0; j < PList_size( gFaces ); j++ ){
                pGFace gFace = (pGFace) PList_item( gFaces , j );
                int  currentTag = GEN_tag( gFace ); // candidate face's tag (debugging only)
                if( gFace != (pGFace) vConG ) { // exclude root classified face
                  GF_parRange ( gFace, 0, &pLow, &pHigh);
                  parFace[0] = ( pLow + pHigh ) * 0.5;
                  GF_parRange ( gFace, 1, &pLow, &pHigh);
                  parFace[1] = ( pLow + pHigh ) * 0.5;
                  echeck = GF_point( (pGFace) gFace , parFace , coordGVOther );
                  dx = coordGVOther[0] - coordGVSelf[0];
                  dy = coordGVOther[1] - coordGVSelf[1];
                  if( dx*dx + dy*dy < FdisTol) {
                    foundETag = GEN_tag( gRegion ); // found Extruded Tag
                    foundEETag = GEN_tag( gFace );   // found Extrusion End Tag
                  }
                }
              }
              PList_delete( gFaces );
              break; 
          } 
          assert(foundEETag != 0);
          count2D++;
          int* vtxData = new int[1];
          vtxData[0] = count2D;
          EN_attachDataPtr((pEntity)vrts[i],myFather,(void*)vtxData);
          V_coord(vrts[i],coordNewPt[i]);

          fprintf ( fcr, "%.15E %.15E %d %d %d %d  \n", coordNewPt[i][0],coordNewPt[i][1], vClassDim, foundESTag, foundETag, foundEETag );
        }
      }

      double coordFather[nvert][3];
      int fatherIds[4]; //store the ids of the fathers (vertices) on the root face
      for(i=0; i< nvert ; i++) {
        int* fatherIdPtr;
        const int exists = EN_getDataPtr((pEntity)vrts[i],myFather,(void**)&fatherIdPtr);
        assert(exists);
        fatherIds[i] = fatherIdPtr[0];
        V_coord(vrts[i],coordFather[i]);
        fprintf ( fcn, "%d ", fatherIds[i]);
      }
      fprintf ( fcn, "\n");

      dir=0;  // 1 fails
      // get the upward adjacent region srcRgn
      region = F_region(face, dir );  // 0 is the negative normal which I assume for a face on the boundary in is interior.
      if(region==NULL) { // try other dir
        dir=1;  // 1 fails
        region = F_region(face, dir );  // 0 is the negative normal which I assume for a face on the boundary in is interior.
      }

      regions=PList_new();
      faces=PList_new();
      err = Extrusion_3DRegionsAndLayerFaces(region, regions, faces, 1);
      PList_delete(regions); // not used so delete
      if(err!=1 && !PCU_Comm_Self())
        fprintf(stderr, "Extrusion_3DRegionsAndLayerFaces returned %d for err \n", err);

      // for each face in the returned list of faces
      iter=0;
      pFace sonFace;
      int iface=0;
      dir=0;
      while( (sonFace = (pFace)PList_next(faces, &iter)) ) { //loop over plist of vertices
        if(iface !=0) {  // root face is in the stack but we already took care of it above
          // get the downward adjacent vertices of face - they will be in the same order as the srcFace ids
          listVn= F_vertices(sonFace, dir);
          void *iter2=0; // Must initialize to 0
          i=0;
          int my2Dfath;
          pVertex  sonVtx;
          double dist, dx, dy, distMin;
          double coordSon[3];
          int iMin;
          while( (sonVtx = (pVertex)PList_next(listVn, &iter2)) ) { //loop over plist of vertices
            V_coord(sonVtx,coordSon);
            distMin=1.0e7;
            for(i=0; i< nvert; i++){
              dx=coordSon[0]-coordFather[i][0];
              dy=coordSon[1]-coordFather[i][1];
              dist=dx*dx+dy*dy;
              if(dist < distMin) {
                iMin=i;
                distMin=dist;
              }
            }
            my2Dfath=fatherIds[iMin];
            int* vtxData = new int[1];
            vtxData[0] = my2Dfath;
            EN_attachDataPtr((pEntity)sonVtx,myFather,(void*)vtxData);
          }
          PList_delete(listVn);
        }
        iface++;
      }
      PList_delete(faces);
    } //end root face iterator
  }
  apf::MeshSIM* cake = reinterpret_cast<apf::MeshSIM*>(simApfMesh);
  cake->createIntTag("fathers2D", myFather, 1);
}

int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  PCU_Comm_Init();
  lion_set_verbosity(1);
  MS_init();
  SimAdvMeshing_start(); //for fancy BL/extrusion queries
  SimModel_start();
  Sim_readLicenseFile(NULL);
  SimPartitionedMesh_start(&argc,&argv);

  getConfig(argc, argv);
  if( should_log )
    Sim_logOn("convert.sim.log");

  if (should_attach_order && should_fix_pyramids) {
    if (!PCU_Comm_Self())
      std::cout << "disabling pyramid fix because --attach-order was given\n";
    should_fix_pyramids = false;
  }

  gmi_sim_start();
  gmi_register_sim();
  pProgress progress = Progress_new();
  Progress_setDefaultCallback(progress);

  gmi_model* mdl;
  if( gmi_native_path )
    mdl = gmi_sim_load(gmi_native_path,gmi_path);
  else
    mdl = gmi_load(gmi_path);

  pGModel simModel = gmi_export_sim(mdl);
/*
  pParasolidNativeModel nModel = ParasolidNM_createFromFile(gmi_native_path,0);
  pGModel    Amodel = GAM_createFromNativeModel(nModel,progress); 
*/
/* leaving this litle model tester in for a while
  pGModel    Amodel = simModel;
  double coordGVOther[3];
  pGVertex gvertex;
  GVIter gvIter=GM_vertexIter(Amodel);
  while ( (gvertex=GVIter_next(gvIter))) {
    int id = GEN_tag(gvertex);
    GV_point( gvertex , coordGVOther );
    cout<< id << " tag and coords " << coordGVOther[0] << " "  << coordGVOther[1]<< " " << coordGVOther[2] <<endl; 
  }
  GVIter_delete(gvIter);
  assert(false);
*/

  double t0 = PCU_Time();
  pParMesh sim_mesh = PM_load(sms_path, simModel, progress);
  double t1 = PCU_Time();
  if(!PCU_Comm_Self())
    fprintf(stderr, "read and created the simmetrix mesh in %f seconds\n", t1-t0);

  apf::Mesh* simApfMesh = apf::createMesh(sim_mesh);

  addFathersTag(simModel, sim_mesh, simApfMesh, extruRootPath);

  double t2 = PCU_Time();
  if(!PCU_Comm_Self())
    fprintf(stderr, "created the apf_sim mesh in %f seconds\n", t2-t1);
  if (should_attach_order) attachOrder(simApfMesh);

  apf::Mesh2* mesh = apf::createMdsMesh(mdl, simApfMesh);
  double t3 = PCU_Time();
  if(!PCU_Comm_Self())
    fprintf(stderr, "created the apf_mds mesh in %f seconds\n", t3-t2);

  apf::printStats(mesh);
  apf::destroyMesh(simApfMesh);
  M_release(sim_mesh);
  fixMatches(mesh);
  if (should_fix_pyramids) fixPyramids(mesh);
  mesh->verify();
  mesh->writeNative(smb_path);

  mesh->destroyNative();
  apf::destroyMesh(mesh);

  Progress_delete(progress);
  gmi_sim_stop();
  SimPartitionedMesh_stop();
  Sim_unregisterAllKeys();
  SimModel_stop();
  MS_exit();
  if( should_log )
    Sim_logOff();
  PCU_Comm_Free();
  MPI_Finalize();
}
