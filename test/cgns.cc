/* Author: Dr Andrew Parker (2019) - FGE Ltd */

#include <apf.h>
#include <gmi_null.h>
#include <gmi_mesh.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <PCU.h>
#include <lionPrint.h>
#include <cstdlib>
//
#include <iostream>
//
#include <parma.h>
#include <apfShape.h>
#include <pumi.h>
#include <apfZoltan.h>

// https://gist.github.com/bgranzow/98087114166956646da684ed98acab02
apf::MeshTag *create_int_tag(const std::string &name, apf::Mesh *m, int dim)
{
  apf::MeshTag *tag = m->createIntTag(name.c_str(), 1); // 1 is size of tag
  apf::MeshEntity *elem = nullptr;
  apf::MeshIterator *it = m->begin(dim);
  int vals[1];
  vals[0] = PCU_Comm_Self();
  while ((elem = m->iterate(it)))
    m->setIntTag(elem, tag, vals);
  m->end(it);
  return tag;
}

//https://github.com/CEED/PUMI/blob/master/ma/maDBG.cc
apf::Field *convert_tag_doubleField(const std::string &name, apf::Mesh *m, apf::MeshTag *t, int dim)
{
  apf::MeshEntity *elem = nullptr;
  apf::MeshIterator *it = m->begin(dim);
  apf::Field *f = nullptr;
  const auto fieldName = name.c_str();
  if (dim == 0)
    f = apf::createFieldOn(m, fieldName, apf::SCALAR);
  else
    f = apf::createField(m, fieldName, apf::SCALAR, apf::getConstant(dim));

  int vals[1];
  while ((elem = m->iterate(it)))
  {
    m->getIntTag(elem, t, vals);
    double dval[1];
    dval[0] = vals[0];
    apf::setComponents(f, elem, 0, dval);
  }
  m->end(it);
  return f;
}

void balance(const std::string &prefix, const apf::ZoltanMethod &method, apf::Mesh2 *m)
{
  const auto dim = m->getDimension();
  convert_tag_doubleField("procID_prebalance", m, create_int_tag("procID_prebalance", m, dim), dim);

  apf::MeshTag *weights = Parma_WeighByMemory(m);
  apf::Balancer *balancer = makeZoltanBalancer(m, method, apf::REPARTITION);
  balancer->balance(weights, 1.10);
  delete balancer;
  apf::removeTagFromDimension(m, weights, m->getDimension());
  Parma_PrintPtnStats(m, "");
  m->destroyTag(weights);

  const std::string name = prefix + "_balance_" + std::to_string(PCU_Comm_Peers()) + "procs";
  apf::writeVtkFiles(name.c_str(), m);
}

//https://github.com/SCOREC/core/blob/4b854ae996cb261a22f6ee6b704569b78866004c/test/reorder.cc
void simpleReorder(const std::string &prefix, apf::Mesh2 *m)
{
  {
    apf::GlobalNumbering *gn = nullptr;
    gn = apf::makeGlobal(apf::numberOwnedNodes(m, "vertex Indices_prereorder"));
    apf::synchronize(gn);
  }
  // no synchronize call
  // https://github.com/SNLComputation/Albany/blob/master/src/disc/pumi/Albany_APFDiscretization.cpp @ various place throughout file
  // https://github.com/SCOREC/core/issues/249
  {
    apf::GlobalNumbering *gn = nullptr;
    gn = apf::makeGlobal(apf::numberElements(m, "element Indices_prereorder"));
  }

  //apf::MeshTag *order = Parma_BfsReorder(m);
  //apf::reorderMdsMesh(m, order);
  apf::reorderMdsMesh(m);

  {
    apf::GlobalNumbering *gn = nullptr;
    gn = apf::makeGlobal(apf::numberOwnedNodes(m, "vertex Indices_postreorder"));
    apf::synchronize(gn);
  }
  // no synchronize call
  // https://github.com/SNLComputation/Albany/blob/master/src/disc/pumi/Albany_APFDiscretization.cpp @ various place throughout file
  // https://github.com/SCOREC/core/issues/249
  {
    apf::GlobalNumbering *gn = nullptr;
    gn = apf::makeGlobal(apf::numberElements(m, "element Indices_postreorder"));
  }

  const std::string name = prefix + "_reorder_" + std::to_string(PCU_Comm_Peers()) + "procs";
  apf::writeVtkFiles(name.c_str(), m);
}

pMesh toPumi(const std::string &prefix, gmi_model *g, apf::Mesh2 *mesh)
{
  //create the pumi instance
  pumi::instance()->model = new gModel(g);
  pMesh pm = pumi_mesh_load(mesh);
  std::cout << pm << std::endl;
  pumi_mesh_verify(pm);
  const std::string name = prefix + "_toPUMI_" + std::to_string(PCU_Comm_Peers()) + "procs";
  pumi_mesh_write(pm, name.c_str(), "vtk");
  return pm;
}

void additional(const std::string &prefix, gmi_model *g, apf::Mesh2 *mesh)
{
  // seems essential to make pm first before calling balance or reorder...
  auto pm = toPumi(prefix, g, mesh);
  balance(prefix, apf::RCB, pm);
  simpleReorder(prefix, pm);

  {
    //create an element field
    const int mdim = pumi_mesh_getDim(pm);
    pShape s = pumi_shape_getConstant(mdim);
    const int dofPerElm = 1;
    pField f = pumi_field_create(pm, "elmField", dofPerElm, PUMI_PACKED, s);

    pMeshIter it = pm->begin(mdim);
    pMeshEnt e;
    double v = 0;
    while ((e = pm->iterate(it)))
      pumi_node_setField(f, e, 0, &v);
    pm->end(it);

    const int ghost = mdim;
    const int bridge = ghost - 1;
    const int numLayers = 1;
    const int ghostAcrossCopies = 1;
    pumi_ghost_createLayer(pm, bridge, ghost, numLayers, ghostAcrossCopies);

    it = pm->begin(mdim);
    v = 1;
    while ((e = pm->iterate(it)))
    {
      if (!pumi_ment_isGhost(e))
        pumi_node_setField(f, e, 0, &v);
    }
    pm->end(it);

    const std::string name = prefix + "_toPUMI_" + std::to_string(PCU_Comm_Peers()) + "procs";
    // only the owned elements will have a elmField value of 1
    pumi_mesh_write(pm, (name + "_beforeSync").c_str(), "vtk");

    pumi_field_synchronize(f);

    // owned and ghosted elements will have a elmField value of 1
    pumi_mesh_write(pm, (name + "_afterSync").c_str(), "vtk");

    // clean-up
    pumi_field_delete(f);
    pumi_ghost_delete(pm);
    pumi_mesh_delete(pm);
  }
}

int main(int argc, char **argv)
{
#ifdef HAVE_CGNS
  MPI_Init(&argc, &argv);
  PCU_Comm_Init();
  lion_set_verbosity(1);
  bool additionalTests = false;
  if (argc < 3)
  {
    if (!PCU_Comm_Self())
      printf("Usage: %s <in .cgns> <out .smb>\n", argv[0]);
    MPI_Finalize();
    exit(EXIT_FAILURE);
    return -1;
  }
  else if (argc == 4)
  {
    if (std::string(argv[3]) == "additional")
      additionalTests = true;
    else
    {
      if (!PCU_Comm_Self())
        printf("Usage: %s <in .cgns> <out .smb> additional\n", argv[0]);
      MPI_Finalize();
      exit(EXIT_FAILURE);
      return -1;
    }
  }
  else if (argc > 4)
  {
    if (!PCU_Comm_Self())
      printf("Usage: %s <in .cgns> <out .smb>\n", argv[0]);
    MPI_Finalize();
    exit(EXIT_FAILURE);
    return -1;
  }

  gmi_register_null();
  gmi_register_mesh();
  gmi_model *g = gmi_load(".null");
  apf::CGNSBCMap cgnsBCMap;
  apf::Mesh2 *m = apf::loadMdsFromCGNS(g, argv[1], cgnsBCMap);
  m->verify();
  //
  m->writeNative(argv[2]);
  // so we can see the result
  const std::string path = argv[1];
  std::size_t found = path.find_last_of("/\\");

#ifndef NDEBUG // debug settings, cmake double negative....
  const auto prefix = path.substr(found + 1) + "_debug";
#else // optimised setting
  const auto prefix = path.substr(found + 1) + "_release";
#endif

  const auto name = prefix + "_" + std::to_string(PCU_Comm_Peers()) + "procs" + std::string("_toVTK");
  std::cout << path << " " << found << " " << name << std::endl;
  apf::writeVtkFiles(name.c_str(), m);

  // main purpose is to call additional tests through the test harness testing.cmake
  if (additionalTests)
    additional(prefix, g, m);

  if (!additionalTests)
  {
    m->destroyNative();
    apf::destroyMesh(m);
  }
  PCU_Comm_Free();
  MPI_Finalize();
  return 0;
#else
  PCU_ALWAYS_ASSERT_VERBOSE(true == false,
                            "Build with ENABLE_CGNS to allow this functionality.");
  exit(EXIT_FAILURE);
  return -1;
#endif
}
