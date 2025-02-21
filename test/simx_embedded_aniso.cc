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

namespace {
  class Args {
  public:
    Args(int argc, char* argv[]);

    #define ARGS_GETTER(name, type) type name(void) const noexcept { return name##_; }
    ARGS_GETTER(help, bool)
    ARGS_GETTER(before, const std::string&)
    ARGS_GETTER(after, const std::string&)
    ARGS_GETTER(input_mesh, const std::string&)
    ARGS_GETTER(input_model, const std::string&)
    ARGS_GETTER(input_nmodel, const std::string&)
    ARGS_GETTER(output_mesh, const std::string&)
    ARGS_GETTER(isotropic, bool)
    #undef ARGS_GETTER

    /** @brief Check for argument parse errors. */
    operator bool() const { return !error_flag_; }

    void print_usage(std::ostream& str) const;
    void print_help(std::ostream& str) const;

  private:
    bool isotropic_{false}, error_flag_{false}, help_{false};
    int verbosity_{0};
    std::string argv0, input_mesh_, input_model_, input_nmodel_, output_mesh_;
    std::string before_, after_;
  }; // class Args
} // namespace

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
  MPI_Init(&argc, &argv);
  PCU_Comm_Init();
  Args args(argc, argv);
  if (args.help()) {
    args.print_help(std::cout);
    PCU_Comm_Free();
    MPI_Finalize();
    return !args ? 1 : 0;
  } else if (!args) {
    args.print_usage(std::cerr);
    PCU_Comm_Free();
    MPI_Finalize();
    return 1;
  }

  cout << endl;
  cout <<" Reading ... " << endl;
  cout <<"  Model from file : " << args.input_model() << endl;
  cout <<"  Mesh from file  : " << args.input_mesh() << endl;

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
  // use pParMesh for compatibility with apf::createMesh
  pParMesh mesh = 0;

  // Load geometry 
  int PARASOLID_TEXT_FORMAT = 0; // 0 for .xmt_txt o .x_t files
  nmodel = !args.input_nmodel().empty() ? ParasolidNM_createFromFile(args.input_nmodel().data(), PARASOLID_TEXT_FORMAT) : NULL;
  if (nmodel && NM_isAssemblyModel(nmodel)) {
    std::cerr << " Parasolid assembly model detected ... cannot handle" << std::endl;
    PCU_Comm_Free();
    MPI_Finalize();
    return 1;
  }
  if (!nmodel) cout<<"  Model is assumed to be discrete"<<endl;
  cout<<endl;

  model = GM_load(args.input_model().data(), nmodel, NULL);
  if(!model) {
    cerr << " ERROR : Didn't load a model, check that code was compiled with the correct MODELER specified and that the model file has an extension such as .xmt_txt, or .XMT_TXT, or .x_t" << endl << endl;
    PCU_Comm_Free();
    MPI_Finalize();
    return 1;
  }
  mesh = PM_load(args.input_mesh().data(), model, NULL);
  
  cout<<endl;
  cout<<" Reading model and mesh done ..."<<endl;
  cout<<endl;

  // Setup and run size field
  ma::Mesh* mesh_ref = apf::createMesh(mesh);
  EmbeddedShockFunction sf(mesh_ref, args.isotropic());

  if (!args.before().empty()) {
    std::cout << " writing size field to before vtk file" << std::endl;

    apf::Field *frameField = nullptr, *scaleField = nullptr;
    frameField = apf::createFieldOn(mesh_ref, "adapt_frames", apf::MATRIX);
    scaleField = apf::createFieldOn(mesh_ref, "adapt_scales", apf::VECTOR);
    ma::Iterator* it = mesh_ref->begin(0);
    for (ma::Entity *v = mesh_ref->iterate(it); v;
      v = mesh_ref->iterate(it)) {
      ma::Vector scale;
      ma::Matrix frame;
      sf.getValue(v, frame, scale);
      apf::setVector(scaleField, v, 0, scale);
      apf::setMatrix(frameField, v, 0, frame);
    }
    mesh_ref->end(it);

    apf::writeVtkFiles(args.before().data(), mesh_ref);
  }
  
  pACase mesh_case = MS_newMeshCase(model);
  MS_setAnisoSizeAttFunc(mesh_case, "anisoUDF", anisoUDF, &sf);
  MS_setAnisoMeshSize(mesh_case, GM_domain(model), MS_userDefinedType | 1, 0, "anisoUDF");
  pMSAdapt adapter = MSA_createFromCase(mesh_case, mesh);
  MSA_adapt(adapter, NULL);

  // Write adapted mesh
  if (!args.output_mesh().empty()) {
    pMesh mesh_write = PM_mesh(mesh,0); // no need to free this according to PM_mesh documentation? 
    cout<<"Adapted mesh statistics: Num. elements: "<<M_numRegions(mesh_write)<<", num. vertices: "<<M_numVertices(mesh_write)<<endl;
    cout<<" start writing adapted mesh: "<<endl;
    M_write(mesh_write,args.output_mesh().data(),0,NULL);
    cout<<" done writing adapted mesh: "<<endl;
  }
  
  if (!args.after().empty()) {
    ma::Mesh* mesh_adapted = apf::createMesh(mesh);
    std::cout << " writing adapted mesh to vtk file" << std::endl;
    apf::writeVtkFiles(args.after().data(), mesh_adapted);
    //apf::destroyMesh(mesh_adapted);
  }

  // Release everything
  //apf::destroyMesh(mesh_ref);
  MS_deleteMeshCase(mesh_case);
  MSA_delete(adapter);
  //M_release(mesh);
  M_release(mesh); 
  GM_release(model);
  NM_release(nmodel);

  // Stop everything
  SimAdvMeshing_stop();  

  SimParasolid_stop(1);
  SimDiscrete_stop(0);

  MS_exit();
  Sim_logOff();

  Sim_unregisterAllKeys();

  PCU_Comm_Free();
  MPI_Finalize();

  return 0;
}

namespace {
  Args::Args(int argc, char* argv[]) {
    argv0 = argv[0];
    int c;
    int given[256] = {0};
    const char* required = "M";
    while ((c = getopt(argc, argv, ":A:B:hM:G:o:vI")) != -1) {
      ++given[c];
      switch (c) {
      case 'A':
        after_ = optarg;
        break;
      case 'B':
        before_ = optarg;
        break;
      case 'h':
        help_ = true;
        break;
      case 'M':
        input_model_ = optarg;
        break;
      case 'G':
        input_nmodel_ = optarg;
        break;
      case 'o':
        output_mesh_ = optarg;
        break;
      case 'v':
        ++verbosity_;
        break;
      case 'I':
        isotropic_ = true;
        break;
      case ':':
        std::cerr << "ERROR: Option -" << char(optopt) << " requires an "
          "argument." << std::endl;
        error_flag_ = true;
        break;
      case '?':
        std::cerr << "ERROR: Unrecognized option: -" << char(optopt) << std::endl;
        error_flag_ = true;
        break;
      }
    }
    for (const char* r = required; *r != '\0'; ++r) {
      if (!given[int(*r)]) {
        std::cerr << "ERROR: Flag -" << *r << " is required." << std::endl;
        error_flag_ = true;
      }
    }
    if (optind < argc) {
      input_mesh_ = argv[optind];
    } else {
      std::cerr << "ERROR: INPUT.sms is required." << std::endl;
      error_flag_ = true;
    }
  } 

  void Args::print_usage(std::ostream& str) const {
    str << "USAGE: " << argv0 << " [-hvI] [-B before.vtk] [-A after.vtk] "
      "[-o OUTPUT.sms] [-M MODEL.smd] [-G MODEL_nat.x_t] INITIAL.sms"
      << std::endl;
  }

  void Args::print_help(std::ostream& str) const {
    print_usage(str);
    str << "simx_aniso adapts a simmetrix mesh with an embedded shock surface.\n";
    str << "OPTIONS:\n"
    "-M MODEL.smd       Set the simmetrix model file (required).\n"
    "-G MODEL_nat.x_t   Set a _nat file, model is discrete when not set.\n"
    "-o OUTPUT.sms      Write final mesh to OUTPUT.sms.\n"
    "-A after.vtk       Write adapted mesh to after.vtk.\n"
    "-B before.vtk      Write initial mesh with size field to before.vtk.\n"
    "-h                 Display this help menu.\n"
    "-I                 Run completely isotropic adaptation for testing purposes.\n"
    "-v                 Increase verbosity. \n";
  }
} // namespace