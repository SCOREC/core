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

//https://github.com/SCOREC/core/blob/4b854ae996cb261a22f6ee6b704569b78866004c/test/repartition.cc
// balance appears to delete the vertex and element numbering
void balance(apf::Mesh2 *m)
{
  bool fineStats = false;                         // set to true for per part stats
  Parma_PrintPtnStats(m, "preRefine", fineStats); //FIXME

  apf::MeshTag *weights = m->createDoubleTag("zoltan_weight", 1);
  {
    apf::MeshIterator *it = m->begin(m->getDimension());
    apf::MeshEntity *e;
    double value = 1.0;
    while ((e = m->iterate(it)))
      m->setDoubleTag(e, weights, &value);
    m->end(it);
  }
  apf::Balancer *b = apf::makeZoltanBalancer(m, apf::RCB, apf::REPARTITION, false);
  b->balance(weights, 1.1);
  m->destroyTag(weights);
  delete b;

  Parma_PrintPtnStats(m, "");
}

//https://github.com/SCOREC/core/blob/4b854ae996cb261a22f6ee6b704569b78866004c/test/reorder.cc
void reorder(apf::Mesh2 *m)
{
  //apf::MeshTag *order = Parma_BfsReorder(m);
  //apf::reorderMdsMesh(m, order);
  apf::reorderMdsMesh(m);
}

// https://gist.github.com/bgranzow/98087114166956646da684ed98acab02
apf::MeshTag *create_int_tag(apf::Mesh *m, int dim)
{
  apf::MeshTag *tag = m->createIntTag("my_tag", 1); // 1 is size of tag
  apf::MeshEntity *elem;
  apf::MeshIterator *it = m->begin(dim);
  int vals[1];
  vals[0] = PCU_Comm_Self();
  while ((elem = m->iterate(it)))
    m->setIntTag(elem, tag, vals);
  m->end(it);
  return tag;
}

//https://github.com/CEED/PUMI/blob/master/ma/maDBG.cc
apf::Field *convert_tag_doubleField(apf::Mesh *m, apf::MeshTag *t, int dim)
{
  apf::MeshEntity *elem;
  apf::MeshIterator *it = m->begin(dim);
  apf::Field *f = nullptr;
  const auto fieldName = "my_field";
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

void additional(gmi_model *g, apf::Mesh2 *mesh)
{
  std::cout << mesh << std::endl;
  {
    balance(mesh);
    const std::string name = "output_balance_" + std::to_string(PCU_Comm_Peers()) + "procs";
    apf::writeVtkFiles(name.c_str(), mesh);
  }
  {
    reorder(mesh);
    const auto dim = mesh->getDimension();
    convert_tag_doubleField(mesh, create_int_tag(mesh, dim), dim);

    {
      apf::GlobalNumbering *gn = nullptr;
      gn = apf::makeGlobal(apf::numberOwnedNodes(mesh, "vertex Indices_postreorder"));
      apf::synchronize(gn);
    }
    // no synchronize call
    // https://github.com/SNLComputation/Albany/blob/master/src/disc/pumi/Albany_APFDiscretization.cpp @ various place throughout file
    // https://github.com/SCOREC/core/issues/249
    {
      apf::GlobalNumbering *gn = nullptr;
      gn = apf::makeGlobal(apf::numberElements(mesh, "element Indices_postreorder"));
    }
    const std::string name = "output_balance_reorder_" + std::to_string(PCU_Comm_Peers()) + "procs";
    apf::writeVtkFiles(name.c_str(), mesh);
  }

  {
    //create the pumi instance
    pumi::instance()->model = new gModel(g);
    pMesh pm = pumi_mesh_load(mesh);
    std::cout << pm << std::endl;
    pumi_mesh_verify(pm);

    // //create an element field
    // const int mdim = pumi_mesh_getDim(pm);
    // pShape s = pumi_shape_getConstant(mdim);
    // const int dofPerElm = 1;
    // pField f = pumi_field_create(pm, "elmField", dofPerElm, PUMI_PACKED, s);

    // pMeshIter it = pm->begin(mdim);
    // pMeshEnt e;
    // double v = 0;
    // while ((e = pm->iterate(it)))
    //   pumi_node_setField(f, e, 0, &v);
    // pm->end(it);

    // const int ghost = mdim;
    // const int bridge = ghost - 1;
    // const int numLayers = 1;
    // const int ghostAcrossCopies = 1;
    // pumi_ghost_createLayer(pm, bridge, ghost, numLayers, ghostAcrossCopies);

    // it = pm->begin(mdim);
    // v = 1;
    // while ((e = pm->iterate(it)))
    // {
    //   if (!pumi_ment_isGhost(e))
    //     pumi_node_setField(f, e, 0, &v);
    // }
    // pm->end(it);

    // // only the owned elements will have a elmField value of 1
    // pumi_mesh_write(pm, "beforeSync", "vtk");

    // pumi_field_synchronize(f);

    // // owned and ghosted elements will have a elmField value of 1
    // pumi_mesh_write(pm, "afterSync", "vtk");

    // // clean-up
    // pumi_field_delete(f);
    // pumi_ghost_delete(pm);
    // pumi_mesh_delete(pm);
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
  apf::Mesh2 *m = apf::loadMdsFromCGNS(argv[1]);
  m->verify();
  //
  m->writeNative(argv[2]);
  // so we can see the result
  const std::string path = argv[1];
  std::size_t found = path.find_last_of("/\\");
  const auto name = path.substr(found + 1) + std::string("_toVTK");
  std::cout << path << " " << found << " " << name << std::endl;
  apf::writeVtkFiles(name.c_str(), m);

  // main purpose is to call additional tests through the test harness testing.cmake
  if (additionalTests)
    additional(g, m);

  m->destroyNative();
  apf::destroyMesh(m);
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
