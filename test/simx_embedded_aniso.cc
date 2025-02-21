#include <stdlib.h>
#include <stdio.h>
#include <strings.h>
#include <string.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <stdexcept>
#include <unistd.h>
#include <getopt.h>

#include <lionPrint.h>
#include <ma.h>
#include <gmi.h>
#include <gmi_sim.h>
#include <apf.h>
#include <apfMesh2.h>
#include <apfZoltan.h>
#include <parma.h>
#include <apfSIM.h>
#include <MeshSim.h>
#include <SimUtil.h>
#include <SimPartitionedMesh.h>

#include "SimParasolidKrnl.h"
#include "MeshSimAdapt.h"
#include "SimDiscrete.h"
#include "SimAdvMeshing.h"
#include "SimMeshTools.h"
#include "embedded_aniso_function.h"

using std::cout;
using std::endl;
using std::cerr;
using std::ofstream;

void anisoUDF(pSizeAttData sadata, void* userdata, double anisosize[3][3]) {
  EmbeddedShockFunction* sf = static_cast<EmbeddedShockFunction*>(userdata);
  double pt[3];
  int haspt = SizeAttData_point(sadata, pt);
  PCU_ALWAYS_ASSERT(haspt);
  
  apf::Vector3 pos;
  pos.fromArray(pt);
  ma::Matrix frame;
  ma::Vector scale;
  sf->getValue(pos, frame, scale);

  // correct size setting (from another working code)
  /*
    anisosize[0][0] = norm[0]*nsize; anisosize[0][1] = norm[1]*nsize; anisosize[0][2] = norm[2]*nsize;
    anisosize[1][0] = tan1[0]*tsize; anisosize[1][1] = tan1[1]*tsize; anisosize[1][2] = tan1[2]*tsize;
    anisosize[2][0] = tan2[0]*tsize; anisosize[2][1] = tan2[1]*tsize; anisosize[2][1] = tan2[2]*tsize; 
  */

  anisosize[0][0] = frame[0][0]*scale[0]; anisosize[0][1] = frame[1][0]*scale[0]; anisosize[0][2] = frame[2][0]*scale[0]; // normal vector
  anisosize[1][0] = frame[0][1]*scale[1]; anisosize[1][1] = frame[1][1]*scale[1]; anisosize[1][2] = frame[2][1]*scale[1]; // tan1
  anisosize[2][0] = frame[0][2]*scale[2]; anisosize[2][1] = frame[1][2]*scale[2]; anisosize[2][1] = frame[2][2]*scale[2]; // tan2

}

int main(int argc, char *argv[])
{
  // Argument parsing
  int displayHelp = 0;
  for(int iArgc=0; iArgc<argc; iArgc++) {
    if(!strcmp(argv[iArgc],"-h")) {
      cout << endl;
      cout << "  HELP requested (by using \"-h\" in arguments) -- " << endl;
      displayHelp = 1;
    }
  }

  if (argc<6 || displayHelp) {
    cout << endl;
    cout << " usage : " << endl;
    cout << "   <executable-name>" << endl;
    cout << "   <model-name.smd>" << endl;
    cout << "   <model-name.x_t or - > ('-' if discrete)" << endl;
    cout << "   <mesh-name.sms>" << endl;
    cout << "   <shock_geometry-model-name.smd or - > ('-' if unused)" << endl;
    cout << "   <model-name.x_t or - > ('-' if unused or if discrete)" << endl;
    cout << endl;
    exit(0);
  }

  char model_file[1024];
  strcpy(model_file,argv[1]);
  char nmodel_file[1024];
  strcpy(nmodel_file,argv[2]);
  char mesh_file[1024];
  strcpy(mesh_file,argv[3]);
  char src_model_file[1024];
  strcpy(src_model_file,argv[4]);
  char src_nmodel_file[1024];
  strcpy(src_nmodel_file,argv[5]);

  cout<<endl;
  cout<<" Reading ... "<<endl;
  cout<<"  Model from file : "<<model_file<<endl;
  cout<<"  Mesh from file  : "<<mesh_file<<endl;

  // Init MeshSim
  Sim_logOn("shock_anisoadapt.log");
  MS_init();
  Sim_readLicenseFile(0);

  SimParasolid_start(1);
  SimDiscrete_start(0);
  SimAdvMeshing_start();

  // Geometry fields
  pNativeModel nmodel = 0;
  pGModel model = 0;
  //pMesh mesh;
  pParMesh mesh;

  // Shock source geometry
  pNativeModel src_nmodel = 0;
  pGModel src_model = 0;

  // Load geometry 
  int PARASOLID_TEXT_FORMAT = 0; // 0 for .xmt_txt o .x_t files
  nmodel = strcmp(nmodel_file, "-") ? ParasolidNM_createFromFile(nmodel_file, PARASOLID_TEXT_FORMAT) : NULL;
  if (nmodel && NM_isAssemblyModel(nmodel)) {
    std::cerr << " Parasolid assembly model detected ... cannot handle" << std::endl;
    exit(0);
  }
  if (!nmodel) cout<<"  Model is assumed to be discrete"<<endl;
  cout<<endl;

  model = GM_load(model_file, nmodel, NULL);
  if(!model) {
    cerr << " ERROR : Didn't load a model, check that code was compiled with the correct MODELER specified and that the model file has an extension such as .xmt_txt, or .XMT_TXT, or .x_t" << endl << endl;
    return 1;
  }
  //mesh = M_load(mesh_file, model, NULL);
  mesh = PM_load(mesh_file, model, NULL); //pParMesh version for apf::createMesh

  // Load source geometry
  if (strcmp(src_model_file, "-")) {
    cout<<" Reading shock source geometry ... "<<endl;
    cout<<"  Model from file : "<<src_model_file<<endl;

    src_nmodel = strcmp(src_nmodel_file, "-") ? ParasolidNM_createFromFile(src_nmodel_file, PARASOLID_TEXT_FORMAT) : NULL;
    if (src_nmodel && NM_isAssemblyModel(src_nmodel)) {
      std::cerr << " Shock source geometry is parasolid assembly model ... cannot handle" << std::endl;
      exit(0);
    }
    if (!src_nmodel) cout<<"  Source model is assumed to be discrete"<<endl;
    cout<<endl;

    src_model = GM_load(src_model_file, src_nmodel, NULL);
    if(!src_model) {
      cerr << " ERROR : Failed to load shock source geometry" << endl << endl;
      return 1;
    }
  }
  
  cout<<endl;
  cout<<" Reading model and mesh done ..."<<endl;
  cout<<endl;

  // Setup and run size field
  ma::Mesh* mesh_ref = apf::createMesh(mesh);
  EmbeddedShockFunction sf(mesh_ref, {});
  
  pACase mesh_case = MS_newMeshCase(model);
  MS_setAnisoSizeAttFunc(mesh_case, "anisoUDF", anisoUDF, &sf);
  MS_setAnisoMeshSize(mesh_case, GM_domain(model), MS_userDefinedType | 1, 0, "anisoUDF");
  pMSAdapt adapter = MSA_createFromCase(mesh_case, mesh);
  MSA_adapt(adapter, NULL);

  // Write adapted mesh
  pMesh mesh_write = PM_mesh(mesh,0);
  cout<<"Adapted mesh statistics: Num. elements: "<<M_numRegions(mesh_write)<<", num. vertices: "<<M_numVertices(mesh_write)<<endl;
  cout<<" start writing adapted mesh: "<<endl;
  M_write(mesh_write,"shock_anisoadapt_output.sms",0,NULL);
  cout<<" done writing adapted mesh: "<<endl;

  // Release everything
  MS_deleteMeshCase(mesh_case);
  MSA_delete(adapter);
  //M_release(mesh);
  M_release(mesh); 
  GM_release(model);
  NM_release(nmodel);

  GM_release(src_model);
  NM_release(src_nmodel);

  // Stop everything
  SimAdvMeshing_stop();  

  SimParasolid_stop(1);
  SimDiscrete_stop(0);

  MS_exit();
  Sim_logOff();

  Sim_unregisterAllKeys();

  return 0;
}