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
#include <exception>
#include <memory>
#include <stdexcept>

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

namespace {

/** \brief Print nested exceptions. */
void print_exception(const std::exception& e, int level = 0);

/** \brief Return the filename without the extension. */
std::string get_no_extension(const std::string& fname) {
  return fname.substr(0, fname.find_last_of("."));
}

/** \brief Print edge cuts. */
void printEdgeCuts(apf::Mesh2 *m);

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

    if (PCU.Self() == 0)
      std::cout << "STATUS: Loading CRE file : " << args.creFilename();
    gmi_model* model = nullptr;
    // Load CRE file sequentially to avoid NFS errors.
    for (int i = 0; i < PCU.Peers(); ++i) {
      if (PCU.Self() == i) {
        if (PCU.Self() == 0) model = gmi_load(args.creFilename().c_str());
        else model = gmi_cap_load_selective(args.creFilename().c_str(), {});
      }
      PCU.Barrier();
    }

    // Create apf::Mesh2 wrapper.
    if (PCU.Self() == 0)
      std::cout << "STATUS: Create apf::Mesh interface" << std::endl;

    ma::Mesh* apfMesh = nullptr;
    bool original = PCU.Self() == 0;
    auto soloPCU = PCU.Split(PCU.Self(), 0);

    if (original) {
      ma::Mesh* apfCapMesh = apf::createCapMesh(model, soloPCU.get());
      apf::disownCapModel(apfCapMesh);
      if (!apf::loadCapSizingFile(
        apfCapMesh, args.sizingFilename(),
        get_no_extension(args.creFilename()) + ".avm.vmap",
        "adapt_scales", "adapt_frames", true, "Kestrel"
      )) {
        throw std::runtime_error("failed to load sizing file");
      }
      apf::writeVtkFiles("before", apfCapMesh);
      /* CONVERT TO APF MDS MESH AND WRITE VTK */
      std::cout << "STATUS: Convert to APF MDS Mesh" << std::endl;
      double t_start = pcu::Time();
      apfMesh = apf::createMdsMesh(model, apfCapMesh, true);
      double t_end = pcu::Time();
      std::cout << "INFO: Mesh converted in " << t_end - t_start << " seconds."
        << std::endl;
      apf::disownMdsModel(apfMesh);
      apf::destroyMesh(apfCapMesh);
      apfMesh->verify();
      apf::printStats(apfMesh);
      apf::writeVtkFiles("before-mds", apfMesh);
    }

    PCU.Barrier();

    parallelAdapt(apfMesh, model, PCU, args.maxiter());

    if (original) {
      apf::Mesh2* mdsMesh = apfMesh;
      apfMesh = apf::makeEmptyCapMesh(model, "MeshAdapt", soloPCU.get());
      apf::disownCapModel(apfMesh);
      auto t0 = pcu::Time();
      /* WRITE MDS MESH TO CRE */
      apf::convert(mdsMesh, apfMesh);
      auto t1 = pcu::Time();
      std::cout << "Converting MDS to CRE: " << t1 - t0 << " seconds"
        << std::endl;
      apf::destroyMesh(mdsMesh);
      
      /* SAVE FINAL CRE FILE */
      auto creOFileName = get_no_extension(args.creFilename())
        + "_adapted.cre";
      gmi_cap_write(model, creOFileName.c_str());
      apf::destroyMesh(apfMesh);
    }
    gmi_destroy(model);
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

