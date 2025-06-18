/**
 * \file capAdapt.cc
 * \brief Test demonstrating adaptation of Capstone mesh based on a bulk sizing
 * field.
 *
 * Sizing field is loaded from a file and vmap. Mesh is converted to MDS
 * temporarily during adapt then exported back to Capstone and saved as CRE.
 * \author Aiden Woodruff
 * \author Eric Mestreau
 */
#include <iostream>
#include <chrono>
#include <exception>
#include <memory>
#include <stdexcept>
#include <thread>
#include <type_traits>

#include <PCU.h>
#include <apf.h>
#include <apfCAP.h>
#include <apfMDS.h>
#include <apfPartition.h>
#include <gmi.h>
#include <gmi_cap.h>
#include <apfConvert.h>
#include <lionPrint.h>
#include <apfMesh2.h>
#include <apfNumbering.h>
#include <ma.h>
#include <pcu_util.h>
#include <parma.h>
#include <apfMETIS.h>
#include <apfZoltan.h>

#include "CapstoneModule.h"
#include "CreateMG_Framework_Mesh.h"
#include "CreateMG_SizingMetricTool.h"

using namespace CreateMG;

namespace {

/** \brief Print nested exceptions. */
void print_exception(const std::exception& e, int level = 0);

/** \brief Load bulk file of one type into std::vector. */
template <class T>
typename std::enable_if<std::is_trivially_copyable<T>::value, void>::type
load_file(const std::string& fname, std::vector<T>& data) {
  try {
    std::ifstream file(fname, std::ios::in | std::ios::binary);
    T val;
    while (file) {
      file.read(reinterpret_cast<char*>(std::addressof(val)), sizeof(T));
      if (file.gcount() == sizeof(T)) data.push_back(val);
    }
  } catch(...) {
    std::throw_with_nested(
      std::runtime_error("failed to load the file `" + fname + "`")
    );
  }
}

/** \brief Return the filename without the extension. */
std::string get_no_extension(const std::string& fname) {
  return fname.substr(0, fname.find_last_of("."));
}

/** \brief Decompose a 6 component metric tensor into frames and scales. */
void decompose_metric(const v_double& metric, apf::Matrix3x3& frame,
  apf::Vector3& scalar) {
  apf::Matrix3x3 m(metric[0], metric[1], metric[2],
    metric[1], metric[3], metric[4],
    metric[2], metric[4], metric[5]);
  int n = apf::eigen(m, &frame[0], &scalar[0]);
  PCU_DEBUG_ASSERT(n == 3);
}

/**
 * \brief Load sizing metric in mesh vertex iteration order.
 * 
 * \param m Capstone mesh database.
 * \param mapFileName path to vmap file with vertex mappings (for ordering).
 * \param sizFielName path to bulk sizing file.
 * */
std::vector<Metric6> loadSizing(
  MDBI *m, std::string mapFileName, std::string sizFileName
);

/** \brief Load Metric6 vector into APF frame and scale fields. */
void convertMetric6toAnisoFields(
  MDBI *m, const std::vector<Metric6>& sizing6,
  apf::Field* frameField, apf::Field* scaleField
);

/** \brief Print edge cuts. */
void printEdgeCuts(apf::Mesh2 *m);

/** \brief Create an apf::Mesh2* from CapstoneModule and bulk sizing file. */
ma::Mesh* makeApfInterfaceWithSizing(
  CapstoneModule& cs, std::unique_ptr<pcu::PCU>& pcu, bool original,
  const std::string& vmapFile, const std::string& sizingFile, bool smooth
);

/** \brief Migrate all elements back to rank 0. */
void migrateHome(apf::Mesh2* mesh);

/** \brief Take a serial MDS mesh, partition, adapt, then localize. */
void parallelAdapt(ma::Mesh* mesh, gmi_model* mod, pcu::PCU& PCU, int maxiter);

/** \brief Command line argument processing class. */
class Args {
public:
  Args() {}
  Args(int argc, char* argv[]) { parse(argc, argv); }
  void parse(int argc, char* argv[]);
  void print_usage(const char* argv0) const;
  /**
   * \brief Get maximum iterations argument.
   *
   * Pass -1 to trust MeshAdapt or a number <= 10 to specify max iterations.
   */
  int maxiter() const noexcept { return maxiter_; }
  /**
   * \brief Get smoothing argument. If true, request smoothing if supported.
   */
  bool smooth() const noexcept { return smooth_; }

  const std::string& creFilename() const noexcept { return creFilename_; }
  const std::string& sizingFilename() const noexcept {
    return sizingFilename_;
  }
private:
  int maxiter_{-1};
  bool smooth_{false};
  std::string creFilename_;
  std::string sizingFilename_;
};

} // namespace

int main(int argc, char** argv) {
  PCU_Init(&argc, &argv);
  { // PCU object scope
  pcu::PCU PCU;
  Args args;
  try {
    args.parse(argc, argv);
  } catch (std::exception& e) {
    std::cerr << "ERROR: ";
    print_exception(e);
    args.print_usage(argv[0]);
    return 1;
  }

  // Enable PUMI library output.
  lion_set_verbosity(1);

  // Start PUMI geometric model interface.
  gmi_cap_start();
  gmi_register_cap();

  try {
    /* LOAD CAPSTONE MESH */
    if (PCU.Self() == 0)
      std::cout << "STATUS: Init Capstone" << std::endl;
    const std::string analysis("Kestrel");
    CapstoneModule cs("capAdapt",
      "Geometry Database : SMLIB",
      "Mesh Database : Create",
      "Attribution Database : Create"
    );
    MDBI *m = cs.get_mesh();
    PCU_ALWAYS_ASSERT(m);

    if (PCU.Self() == 0)
      std::cout << "STATUS: Loading CRE file : " << args.creFilename()
        << std::endl;
    M_GModel gmodel;
    // Load CRE file sequentially to avoid NFS errors.
    for (int p = 0; p < PCU.Peers(); ++p) {
      if (PCU.Self() == p) {
        gmodel = cs.load_files(v_string(1, args.creFilename()));
      }
      PCU.Barrier();
    }

    // Load mesh.
    M_MModel mmodel;
    std::vector<M_MModel> mmodels;
    MG_API_CALL(m, get_associated_mesh_models(gmodel, mmodels));
    PCU_ALWAYS_ASSERT(mmodels.size() == 1);
    MG_API_CALL(m, set_current_model(mmodels[0]));
    mmodel = mmodels[0];
    cs.set_analysis_name(analysis);

    /* SET THE ADJACENCIES */
    if (PCU.Self() == 0) {
      std::cout << "STATUS: Compute required adjacencies" << std::endl;
      MG_API_CALL(m, set_adjacency_state(REGION2FACE | REGION2EDGE |
        REGION2VERTEX | FACE2EDGE | FACE2VERTEX));
      MG_API_CALL(m, set_reverse_states());
      MG_API_CALL(m, compute_adjacency());
    }

    // Create apf::Mesh2 wrapper.
    if (PCU.Self() == 0)
      std::cout << "STATUS: Create apf::Mesh interface" << std::endl;

    ma::Mesh* apfCapMesh = nullptr, *apfMesh = nullptr;
    bool original = PCU.Self() == 0;
    int parts = PCU.Peers();
    auto soloPCU = PCU.Split(PCU.Self(), 0);

    apfCapMesh = makeApfInterfaceWithSizing(
      cs, soloPCU, original,
      get_no_extension(args.creFilename()) + ".avm.vmap",
      args.sizingFilename(), args.smooth()
    );
    if (original) {
      /* CONVERT TO APF MDS MESH AND WRITE VTK */
      std::cout << "STATUS: Convert to APF MDS Mesh" << std::endl;
      double t_start = pcu::Time();
      apfMesh = apf::createMdsMesh(apfCapMesh->getModel(), apfCapMesh, true);
      double t_end = pcu::Time();
      std::cout << "INFO: Mesh converted in " << t_end - t_start << " seconds."
        << std::endl;
      apfMesh->verify();
      apf::printStats(apfMesh);
      apf::writeVtkFiles("before-mds", apfMesh);
    }

    PCU.Barrier();

    parallelAdapt(apfMesh, apfCapMesh->getModel(), PCU, args.maxiter());

    if (original) {
      /* COPY THE MESH MODEL TO KEEP */
      MG_API_CALL(m, create_associated_model(mmodel, gmodel, "MeshAdapt"));

      /* SETUP NEW MESH ADJACENCIES */
      MG_API_CALL(m, set_adjacency_state(REGION2FACE|REGION2EDGE|
        REGION2VERTEX|FACE2EDGE|FACE2VERTEX));
      MG_API_CALL(m, set_reverse_states());
      MG_API_CALL(m, compute_adjacency());

      // Destroy old fields.
      apf::destroyField(apfCapMesh->findField("adapt_scales"));
      apf::destroyField(apfCapMesh->findField("adapt_frames"));

      auto t0 = pcu::Time();
      /* WRITE MDS MESH TO CRE */
      apf::convert(apfMesh, apfCapMesh);
      auto t1 = pcu::Time();
      std::cout << "Converting MDS to CRE: " << t1 - t0 << " seconds"
        << std::endl;
      
      /* SAVE FINAL CRE FILE */
      cs.save_mesh_file("after_cap.vtk", mmodel);
      auto creOFileName = get_no_extension(args.creFilename())
        + "_adapted.cre";
      cs.save_file(creOFileName,gmodel);
      
      /* PRINT ADAPTED MESH INFO */
      std::string info;
      m->print_info(mmodel, info);
      std::cout << info << std::endl;

      apf::destroyMesh(apfCapMesh);
      apf::destroyMesh(apfMesh);
    }
    PCU.Barrier();
  } catch (const std::exception& e) {
    std::cerr << "FATAL: ";
    print_exception(e);
  } catch(...) {
    std::cerr << "FATAL: unspecified error" << std::endl;
  }

  gmi_cap_stop();
  } // PCU object scope
  PCU_Finalize();
  return 0;
}

namespace {

void print_exception(const std::exception& e, int level) {
  std::cerr << std::string(level * 2, ' ') << e.what() << '\n';
  try {
    std::rethrow_if_nested(e);
  } catch (const std::exception& nestedE) {
    print_exception(nestedE, level + 1);
  } catch (...) {}
}

std::vector<Metric6> loadSizing(
  MDBI *m, std::string mapFileName, std::string sizFileName
) {
  std::vector<std::size_t> vmap;
  std::vector<double>      sizing;
  std::cout << "STATUS: Loading VMap file: " << mapFileName << std::endl;
  load_file<std::size_t>(mapFileName, vmap);
  std::cout << "INFO: VMap number of vertices = " << vmap.size() - 1
    << "    vmap[0] = " << vmap[0]  << std::endl;
  std::cout << "STATUS: Loading Sizing File: " << sizFileName
    << std::endl;
  load_file<double>(sizFileName, sizing);
  std::cout << "INFO: Sizing File number of vertices = "
    << sizing.size()/6 << std::endl;
  std::size_t numVertices;
  MG_API_CALL(m, get_num_topos(Mesh::TOPO_VERTEX, numVertices));
  if (sizing.size() / 6 != numVertices) {
    throw std::runtime_error(
      "sizing field does not match mesh: " + std::to_string(sizing.size() / 6)
      + " / " + std::to_string(numVertices)
    );
  }

  /* COPY SIZEFIELD TO METRIC6 */
  std::vector<Metric6> sizing6(numVertices);
  for (size_t si = 0, i = 1; i < vmap.size(); ++i) {
    size_t tid = vmap[i];
    PCU_DEBUG_ASSERT(tid != 0);
    M_MTopo mtopo = m->get_topo_by_id(Mesh::TOPO_VERTEX, tid);
    for (size_t j = 0; j < 6; ++j, ++si) {
      sizing6[tid - 1][j] = sizing[si];
    }
  }
  return sizing6;
}

void convertMetric6toAnisoFields(
  MDBI *m, const std::vector<Metric6>& sizing6,
  apf::Field* frameField, apf::Field* scaleField
) {
  v_double       sz(6);
  apf::Vector3   scalar;
  apf::Matrix3x3 frame;
  for (std::size_t i=0; i<sizing6.size(); i++) {
    M_MTopo mtopo = m->get_topo_by_id(CreateMG::Mesh::TOPO_VERTEX,i + 1);
    for (std::size_t j=0; j<6; j++) sz[j] = sizing6[i][j];
    
    // field eigen values
    decompose_metric(sz,frame,scalar);
    for (int i=0; i<3; i++) {
      scalar[i] = 1.0/sqrt(scalar[i]);
    }
    // MeshAdapt wants frames on the columns.
    frame = apf::transpose(frame);
    apf::MeshEntity* ent = (apf::MeshEntity*)mtopo.get();
    apf::setVector(scaleField, ent, 0, scalar);
    apf::setMatrix(frameField, ent, 0, frame);
  }
}

ma::Mesh* makeApfInterfaceWithSizing(
  CapstoneModule& cs, std::unique_ptr<pcu::PCU>& pcu, bool original,
  const std::string& vmapFile, const std::string& sizingFile, bool smooth
) {
  // All ranks need apfCapMesh->getModel().
  ma::Mesh* apfCapMesh = apf::createCapMesh(
    cs.get_mesh(), cs.get_geometry(), pcu.get()
  );

  if (original) {
    /* LOAD THE SIZING FILE AND VMAP FILE                        */
    /* CHECK VALIDITY OF THE SIZING FIELD AGAINST THE MESH MODEL */
    std::cout << "STATUS: Loading sizing" << std::endl;
    std::vector<Metric6> sizing6 = loadSizing(
      cs.get_mesh(), vmapFile, sizingFile
    );

    if (smooth) {
      /* SMOOTH SIZING FIELD */
      std::cout << "STATUS: Smoothing sizing field." << std::endl;
      auto smooth_tool = get_sizing_metric_tool(
        cs.get_context(), "CreateSmoothingBase"
      );
      if (smooth_tool == nullptr) {
        throw std::runtime_error("unable to find \"CreateSmoothingBase\"");
      } else {
        smooth_tool->set_context(cs.get_context());
        M_MModel mmodel;
        MG_API_CALL(cs.get_mesh(), get_current_model(mmodel));
        smooth_tool->set_metric(mmodel, "sizing6", sizing6);
        std::vector<Metric6> ometric;
        auto analysis = cs.get_analysis_name();
        smooth_tool->smooth_metric(mmodel, analysis, "sizing6", ometric);
        sizing6 = ometric;
      }
    } else {
      std::cout << "INFO: not smoothing sizing field." << std::endl;
    }

    /* CREATE THE ANISOTROPIC FIELD */
    std::cout << "STATUS: Creating Anisotropic size-field from the file: "
      << sizingFile << std::endl;
    apf::Field* frameField = apf::createFieldOn(
      apfCapMesh, "adapt_frames", apf::MATRIX
    );
    apf::Field* scaleField = apf::createFieldOn(
      apfCapMesh, "adapt_scales", apf::VECTOR
    );
    convertMetric6toAnisoFields(
      cs.get_mesh(), sizing6, frameField, scaleField
    );

    /* WRITE THE BEFORE ADAPT MESH TO VTK USING APF VTK WRITER */
    std::cout << "STATUS: Writing VTK file (before)" << std::endl;
    apf::writeVtkFiles("before", apfCapMesh);
  }
  return apfCapMesh;
}

void migrateHome(apf::Mesh2* mesh) {
  auto t0 = pcu::Time();
  apf::Migration* plan = new apf::Migration(mesh);
  apf::MeshIterator* it = mesh->begin(mesh->getDimension());
  for (apf::MeshEntity* e = mesh->iterate(it); e;
    e = mesh->iterate(it)) {
    plan->send(e, 0);
  }
  mesh->end(it);
  mesh->migrate(plan); // destroys plan
  struct Map0 : apf::Remap {
    virtual int operator()(int) { return 0; }
  } map0;
  apf::remapPartition(mesh, map0);
  auto t1 = pcu::Time();
  if (mesh->getPCU()->Self() == 0)
    std::cout << "INFO: Migrated back to original rank in " << t1 - t0
      << " seconds" << std::endl;
}

void parallelAdapt(ma::Mesh* mesh, gmi_model* model, pcu::PCU& PCU, int maxiter) {
  auto t0 = pcu::Time();
  bool original = false;
  pcu::PCU* oldPCU = nullptr;
  if (mesh) {
    original = true;
    oldPCU = mesh->getPCU();
  }
  int parts = PCU.Peers();
  apf::Migration* plan = nullptr;

  if (original && parts > 1) {
    std::cout << "STATUS: Partitioning the mesh." << std::endl;
    #if defined(PUMI_HAS_ZOLTAN)
    apf::Splitter* splitter = apf::makeZoltanSplitter(
      mesh, apf::GRAPH, apf::PARTITION
    );
    #elif defined(PUMI_HAS_METIS)
    apf::Splitter* splitter = apf::makeMETISsplitter(mesh);
    #else
    apf::Splitter* splitter = Parma_MakeRibSplitter(mesh);
    #endif
    apf::MeshTag* weights = Parma_WeighByMemory(mesh);
    // Split into pieces with 10% imbalance.
    plan = splitter->split(weights, 1.10, parts);
    apf::removeTagFromDimension(mesh, weights,
      mesh->getDimension());
    mesh->destroyTag(weights);
    delete splitter;
    mesh->switchPCU(&PCU);
  }
  PCU.Barrier();
  if (parts > 1) {
    mesh = apf::repeatMdsMesh(
      mesh, model, plan, parts, &PCU
    );
    plan = nullptr; // plan is freed by apf::repeatMdsMesh.
    apf::printStats(mesh);
    apf::writeVtkFiles("before-mds-split", mesh);
  }
  apf::disownMdsModel(mesh);
  printEdgeCuts(mesh);

  apf::Field* mdsScaleField = mesh->findField("adapt_scales");
  apf::Field* mdsFrameField = mesh->findField("adapt_frames");
  PCU_ALWAYS_ASSERT(mdsScaleField);
  PCU_ALWAYS_ASSERT(mdsFrameField);

  /* SETUP AND CALL ADAPT */
  if (original) std::cout << "STATUS: Adaptation setup" << std::endl;
  // adapt setup
  ma::Input* in = ma::makeAdvanced(
    ma::configure(mesh, mdsScaleField, mdsFrameField)
  );
  in->shouldForceAdaptation = true;
  if (maxiter != -1) in->maximumIterations = maxiter;

  if (original) std::cout << "STATUS: Adapting" << std::endl;
  ma::adapt(in);

  /* WRITE THE AFTER ADAPT MESH TO VTK USING APF VTK WRITER */
  if (original) std::cout << "STATUS: Writing VTK file (after)" << std::endl;
  apf::writeVtkFiles("after-mds", mesh);

  if (parts > 1) migrateHome(mesh);
  if (original) {
    auto t1 = pcu::Time();
    std::cout << "INFO: MDS end-to-end time: " << t1 - t0 << " seconds"
      << std::endl;
    mesh->switchPCU(oldPCU);
  } else {
    apf::destroyMesh(mesh);
  }
}

void printEdgeCuts(apf::Mesh2 *m) {
  int ownedFaces = 0, sharedFaces = 0, ownedShared = 0, ct = m->count(2);
  apf::MeshIterator *it = m->begin(2);
  for (apf::MeshEntity *e = m->iterate(it); e; e = m->iterate(it)) {
    if (m->isShared(e)) ++sharedFaces;
    if (m->isOwned(e)) ++ownedFaces;
    if (m->isOwned(e) && m->isShared(e)) ++ownedShared;
  }
  m->end(it);
  int peers = m->getPCU()->Peers();
  std::vector<int> allOwned(peers), allShared(peers), allCt(peers),
    allOwnedShared(peers);
  m->getPCU()->Allgather(&ownedFaces, allOwned.data(), 1);
  m->getPCU()->Allgather(&sharedFaces, allShared.data(), 1);
  m->getPCU()->Allgather(&ownedShared, allOwnedShared.data(), 1);
  m->getPCU()->Allgather(&ct, allCt.data(), 1);
  if (m->getPCU()->Self() == 0) {
    for (int i = 0; i < peers; ++i) {
      std::cout << "Part " << i << ": "
        << allCt[i] << " faces, "
        << allShared[i] << " shared faces, "
        << allOwned[i] << " owned faces, "
        << allOwnedShared[i] << " owned and shared." << std::endl;
    }
  }
}

void Args::parse(int argc, char* argv[]) {
  for (int i = 1; i < argc; ++i) {
    if (argv[i][0] == '-') {
      switch (argv[i][1]) {
      case 'i':
        if (argv[i][2]) maxiter_ = std::stoi(&argv[i][2]);
        else if (++i < argc) maxiter_ = std::stoi(argv[i]);
        else throw std::invalid_argument("missing argument to -i");
        break;
      case 's':
        smooth_ = true;
        break;
      default:
        throw std::invalid_argument("invalid option -" + argv[i][1]);
      }
    } else {
      if (creFilename_.empty()) creFilename_ = argv[i];
      else if (sizingFilename_.empty()) sizingFilename_ = argv[i];
      else {
        throw std::invalid_argument(
          std::string("invalid argument `") + argv[i] + "`"
        );
      }
    }
  }
}

void Args::print_usage(const char* argv0) const {
  std::cout << "USAGE: " << argv0 << " [-i MAXITER] [-s] INPUT.CRE SIZING.DAT"
    << std::endl;
}

} // namsepace

