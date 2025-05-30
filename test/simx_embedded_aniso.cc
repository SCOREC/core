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
#include <PCU.h>

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
    ARGS_GETTER(mesh_adapt, const std::string&)
    ARGS_GETTER(isotropic, bool)
    ARGS_GETTER(planar_opt, int)
    ARGS_GETTER(global_size, double)
    ARGS_GETTER(norm_size, double)
    ARGS_GETTER(ref_radius, double)
    ARGS_GETTER(shock_thickness, double)
    ARGS_GETTER(tip, double)
    ARGS_GETTER(upstream, double)
    ARGS_GETTER(verbosity, bool)
    ARGS_GETTER(mds_maxiter, int)
    ARGS_GETTER(refine_only, bool)
    #undef ARGS_GETTER

    /** @brief Check for argument parse errors. */
    operator bool() const { return !error_flag_; }

    void print_usage(std::ostream& str) const;
    void print_help(std::ostream& str) const;

  private:
    bool isotropic_{false}, error_flag_{false}, help_{false}, refine_only_{false};
    int verbosity_{0}, mds_maxiter_{-1}, planar_opt_{-1};
    double global_size_{-1}, norm_size_{-1}, ref_radius_{-1}, shock_thickness_{-1}, tip_{-1}, upstream_{-1};
    std::string argv0, input_mesh_, input_model_, input_nmodel_, output_mesh_;
    std::string before_, after_;
    std::string mesh_adapt_;
  }; // class Args
} // namespace

void anisoUDF(pSizeAttData sadata, void* userdata, double anisosize[3][3]) {
  EmbeddedShockFunction* sf = static_cast<EmbeddedShockFunction*>(userdata);
  double pt[3];
  int haspt = SizeAttData_point(sadata, pt);
  PCU_ALWAYS_ASSERT(haspt);
  
  apf::Vector3 pos; 
  pos.fromArray(pt);
  sf->getSimValue(pos, anisosize);
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

  lion_set_verbosity(args.verbosity());

  cout << endl;
  cout <<" Reading ... " << endl;
  cout <<"  Model from file : " << args.input_model() << endl;
  cout <<"  Mesh from file  : " << args.input_mesh() << endl;

  // Init MeshSim
  Sim_logOn("shock_anisoadapt.log");
  MS_init();
  SimAdvMeshing_start();
  SimModel_start();
  Sim_readLicenseFile(NULL);
  SimPartitionedMesh_start(&argc,&argv);

  gmi_sim_start();
  gmi_register_sim();
  pProgress progress = Progress_new();
  Progress_setDefaultCallback(progress);

  // Simmetrix
  //pNativeModel nmodel = 0;
  pGModel model = 0;
  pParMesh mesh = 0;

  // SCOREC
  gmi_model* mdl_ref;
  ma::Mesh* mesh_ref;

  // Load geometry (closely following convert.cc)
  // todo: maybe some modules not loaded?
  if(!args.input_nmodel().empty()) {
    mdl_ref = gmi_sim_load(args.input_nmodel().data(),args.input_model().data());
  } else {
    mdl_ref = gmi_load(args.input_model().data());
  }
  model = gmi_export_sim(mdl_ref);

  // Load mesh
  mesh = PM_load(args.input_mesh().data(), model, progress);
  mesh_ref = apf::createMesh(mesh);
  
  cout<<endl;
  cout<<" Reading model and mesh done ..."<<endl;
  cout<<"Initial mesh statistics: Num. elements: "<<M_numRegions(PM_mesh(mesh,0))<<", num. vertices: "<<M_numRegions(PM_mesh(mesh,0))<<endl;
  cout<<endl;

  // Setup and run size field
  AnalyticClosestPoint closestPointFunction = doubleConeClosestPointAnalytic;
  if(args.planar_opt() == 0) {
    std::cout << " trying planar shock (option 0)" << std::endl;
    closestPointFunction = planarClosestPointAnalytic;
  } else if (args.planar_opt() == 1) {
    std::cout << " trying angled planar shock (option 1)" << std::endl;
    closestPointFunction = planar30DegTowardsY;
  }
  EmbeddedShockFunction sf(mesh_ref, args.isotropic(), args.global_size(), args.norm_size(), args.shock_thickness(), args.ref_radius(), \
    args.tip(), args.upstream(), closestPointFunction);

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

  // Try SCOREC meshadapt
  if (!args.mesh_adapt().empty()) {

    cout << " trying meshadapt" << endl;

    apf::Mesh2* mesh_mds = apf::createMdsMesh(mdl_ref, mesh_ref);
    EmbeddedShockFunction sf_simx(mesh_mds, args.isotropic(), args.global_size(), args.norm_size(), args.shock_thickness(), args.ref_radius(), \
      args.tip(), args.upstream(), closestPointFunction);
    ma::Input* in = ma::makeAdvanced(ma::configure(mesh_mds, &sf_simx));

    if (args.mds_maxiter() >= 0) {
      cout << " maximum mds adapt iterations set to " << args.mds_maxiter() << endl;
      in->maximumIterations = args.mds_maxiter();
    }
    if(args.refine_only()) {
      std::cout << "Doing refinement only (no coarsening, snapping, or shape fix)." << std::endl;
      in->shouldCoarsen = false;
      in->shouldSnap=false;
      in->shouldFixShape=false;
    }

    double t0 = PCU_Time();
    if (args.verbosity() > 1) {
      ma::adaptVerbose(in, args.verbosity() > 2);
    } else {
      ma::adapt(in);
    }
    cout << " MeshAdapt adapt time (s): " << PCU_Time()-t0 << endl;

    cout<<"Adapted mesh statistics (MDS): Num. elements: "<<mesh_mds->count(3)<<", num. vertices: "<<mesh_mds->count(0)<<endl;
    cout << " writing meshdapt result" << endl;
    apf::writeVtkFiles(args.mesh_adapt().data(), mesh_mds);

    cout << " destroying meshdapt objects" << endl;
    apf::destroyMesh(mesh_mds);
  }

  if (!args.output_mesh().empty() || !args.after().empty()) {
    cout << " trying simmetrix adapt" << endl;

    pACase mesh_case = MS_newMeshCase(model);
    MS_setAnisoSizeAttFunc(mesh_case, "anisoUDF", anisoUDF, &sf);
    MS_setAnisoMeshSize(mesh_case, GM_domain(model), MS_userDefinedType | 1, 0, "anisoUDF");
    pMSAdapt adapter = MSA_createFromCase(mesh_case, mesh);
    MSA_setNewVertexSizeMode(adapter, 1);

    // Try setting size for each vertex before adapts
    /*
    pAllEntProcIter eit = PM_allEntProcIter(mesh, 0);
    pEntity ent;
    pVertex vtx;
    while(ent = AllEntProcIter_next(eit)) {
      vtx = static_cast<pVertex>(ent);
      double anisosize[3][3];

      double pt[3];
      V_coord(vtx, pt);
      apf::Vector3 pos; 
      pos.fromArray(pt);

      sf.getSimValue(pos, anisosize);
      MSA_setAnisoVertexSize(adapter, vtx, anisosize);
    }
    AllEntProcIter_delete(eit);
    */

    double t1 = PCU_Time();
    MSA_adapt(adapter, NULL);
    cout << " Simmetrix adapt time (s): " << PCU_Time()-t1 << endl;

    if (!args.output_mesh().empty()) {
      pMesh mesh_write = PM_mesh(mesh,0); // no need to free this according to PM_mesh documentation? 
      cout<<"Adapted mesh statistics (Simmetrix): Num. elements: "<<M_numRegions(mesh_write)<<", num. vertices: "<<M_numVertices(mesh_write)<<endl;
      cout<<" start writing adapted mesh: "<<endl;
      M_write(mesh_write,args.output_mesh().data(),0,NULL);
      cout<<" done writing adapted mesh: "<<endl;
    }

    if (!args.after().empty()) {
      ma::Mesh* mesh_adapted = apf::createMesh(mesh);
      cout << " writing adapted mesh to vtk file" << endl;
      apf::writeVtkFiles(args.after().data(), mesh_adapted);
    }

    MS_deleteMeshCase(mesh_case);
    MSA_delete(adapter);
  }

  // Release everything
  M_release(mesh); 
  //GM_release(model);

  // Stop everything
  /*
  SimAdvMeshing_stop();  
  SimParasolid_stop(1);
  SimDiscrete_stop(0);
  SimPartitionedMesh_stop();
  SimModel_stop();
  */

  Progress_delete(progress);
  gmi_sim_stop();
  SimPartitionedMesh_stop();
  Sim_unregisterAllKeys();
  SimModel_stop();

  MS_exit();
  Sim_logOff();
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
    while ((c = getopt(argc, argv, ":A:B:hM:G:o:vIPm:i:g:n:r:t:Rp:u:")) != -1) {
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
      case 'm':
        mesh_adapt_ = optarg;
        break;
      case 'v':
        ++verbosity_;
        break;
      case 'I':
        isotropic_ = true;
        break;
      case 'P':
        ++planar_opt_;
        break;
      case 'g':
        global_size_ = std::atof(optarg);
        break;
      case 'n':
        norm_size_ = std::atof(optarg);
        break;
      case 'r':
        ref_radius_ = std::atof(optarg);
        break;
      case 't':
        shock_thickness_ = std::atof(optarg);
        break;
      case 'p':
        tip_ = std::atof(optarg);
        break;
      case 'u':
        upstream_ = std::atof(optarg);
        break;
      case 'i':
        mds_maxiter_ = std::atoi(optarg);
        break;
      case 'R':
        refine_only_ = true;
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
    str << "USAGE: " << argv0 << " [-hvIR] [-B before.vtk] [-A simxadapt.vtk] "
      "[-o SIMX_OUTPUT.sms] [-M MODEL.smd] [-G MODEL_nat.x_t] "
      "[-m meshadapt.vtk] [-i mds maxiter.] "
      "[-g h_global] [-n h_norm] [-r tr_radius] INITIAL.sms"
      << std::endl;
  }

  void Args::print_help(std::ostream& str) const {
    print_usage(str);
    str << "simx_aniso adapts a simmetrix mesh with an embedded shock surface.\n";
    str << "OPTIONS:\n"
    "-M MODEL.smd       Set the simmetrix model file (required).\n"
    "-G MODEL_nat.x_t   Set a _nat file, model is discrete when not set.\n"
    "-o SIMX_OUTPUT.sms Try Simmetrix adapt and write final mesh to OUTPUT.sms.\n"
    "-A simxadapt.vtk   Write Simmetrix adapted mesh to simxadapt.vtk.\n"
    "-B before.vtk      Write initial mesh with size field to before.vtk.\n"
    "-m meshadapt.vtk   Try SCOREC MeshAdapt and write to meshadapt.vtk.\n"
    "-i mds MAXITER     Maximum MDS adapt iterations.\n"
    "-R                 Run refinement only (disable coarsening, snapping and shape fix)\n"
    "-h                 Display this help menu.\n"
    "-I                 Run completely isotropic adaptation for testing purposes.\n"
    "-P                 Planar shock for testing."
    "-g h_global        Override h_global.\n"
    "-n h_norm          Override h_norm.\n"
    "-r tr_radius       Override tip refinement radius.\n"
    "-t thickness       Override thickness.\n"
    "-p h_tip           Override h_tip.\n"
    "-u h_upstream      Override h_upstream.\n"
    "-v                 Increase verbosity. \n";
  }
} // namespace