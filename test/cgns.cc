/* Author: Dr Andrew Parker (2019) - FGE Ltd */

#include <apf.h>
#include <gmi_null.h>
#include <gmi_mesh.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <PCU_C.h>
#include <lionPrint.h>
#include <cstdlib>
//
#include <iostream>
#include <functional>
//
#include <parma.h>
#include <apfShape.h>
#include <pumi.h>
#include <apfZoltan.h>

// https://gist.github.com/bgranzow/98087114166956646da684ed98acab02
apf::MeshTag *create_int_tag(const std::string &name, apf::Mesh2 *m, int dim)
{
  apf::MeshTag *tag = m->createIntTag(name.c_str(), 1); // 1 is size of tag
  apf::MeshEntity *elem = nullptr;
  apf::MeshIterator *it = m->begin(dim);
  int vals[1];
  vals[0] = m->getPCU()->Self();
  while ((elem = m->iterate(it)))
    m->setIntTag(elem, tag, vals);
  m->end(it);
  return tag;
}

//https://github.com/CEED/PUMI/blob/master/ma/maDBG.cc
apf::Field *convert_tag_doubleField(const std::string &name, apf::Mesh2 *m, apf::MeshTag *t, int dim)
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
    double dval;
    dval = vals[0];
    apf::setScalar(f, elem, 0, dval);
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

  const std::string name = prefix + "_balance_" + std::to_string(m->getPCU()->Peers()) + "procs";
  apf::writeVtkFiles(name.c_str(), m);
}

//https://github.com/SCOREC/core/blob/4b854ae996cb261a22f6ee6b704569b78866004c/test/reorder.cc
void simpleReorder(const std::string &prefix, apf::Mesh2 *m)
{
  {
    apf::GlobalNumbering *gn = nullptr;
    gn = apf::makeGlobal(apf::numberOwnedNodes(m, "vert Idx_preOrder"));
    apf::synchronize(gn);
  }
  // no synchronize call
  // https://github.com/SNLComputation/Albany/blob/master/src/disc/pumi/Albany_APFDiscretization.cpp @ various place throughout file
  // https://github.com/SCOREC/core/issues/249
  {
    apf::makeGlobal(apf::numberElements(m, "elem Idx_preOrder"));
  }

  //apf::MeshTag *order = Parma_BfsReorder(m);
  //apf::reorderMdsMesh(m, order);
  apf::reorderMdsMesh(m);

  {
    apf::GlobalNumbering *gn = nullptr;
    gn = apf::makeGlobal(apf::numberOwnedNodes(m, "vert Idx_pstOrder"));
    apf::synchronize(gn);
  }
  // no synchronize call
  // https://github.com/SNLComputation/Albany/blob/master/src/disc/pumi/Albany_APFDiscretization.cpp @ various place throughout file
  // https://github.com/SCOREC/core/issues/249
  {
    apf::makeGlobal(apf::numberElements(m, "elem Idx_pstOrder"));
  }

  const std::string name = prefix + "_reorder_" + std::to_string(m->getPCU()->Peers()) + "procs";
  apf::writeVtkFiles(name.c_str(), m);
}

pMesh toPumi(const std::string &prefix, gmi_model *g, apf::Mesh2 *mesh)
{
  //create the pumi instance
  pumi::instance()->model = new gModel(g);
  pMesh pm = pumi_mesh_load(mesh);
  std::cout << pm << std::endl;
  pumi_mesh_verify(pm);
  const std::string name = prefix + "_toPUMI_" + std::to_string(pm->getPCU()->Peers()) + "procs";
  pumi_mesh_write(pm, name.c_str(), "vtk");
  return pm;
}

auto additional(const std::string &prefix, gmi_model *g, apf::Mesh2 *mesh)
{
  // seems essential to make pm first before calling balance or reorder...
  auto pm = toPumi(prefix, g, mesh);
  balance(prefix, apf::RCB, pm);
  simpleReorder(prefix, pm);

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

  const std::string name = prefix + "_toPUMI_" + std::to_string(mesh->getPCU()->Peers()) + "procs";
  // only the owned elements will have a elmField value of 1
  pumi_mesh_write(pm, (name + "_beforeSync").c_str(), "vtk");

  pumi_field_synchronize(f);

  // owned and ghosted elements will have a elmField value of 1
  pumi_mesh_write(pm, (name + "_afterSync").c_str(), "vtk");

  const auto clean = [&pm, &f]() {
    // clean-up
    pumi_field_delete(f);
    pumi_ghost_delete(pm);
    pumi_mesh_delete(pm);
  };
  return clean;
}

std::string doit(apf::CGNSBCMap &cgnsBCMap, const std::string &argv1, const std::string &argv2, const std::string &post, const bool &additionalTests, const bool &writeCGNS, const std::vector<std::pair<std::string, std::string>> &meshData, pcu::PCU *PCUObj)
{
  gmi_register_null();
  gmi_register_mesh();
  gmi_model *g = gmi_load(".null");
  PCU_t h;
  h.ptr = static_cast<void*>(PCUObj);
  apf::Mesh2 *m = apf::loadMdsFromCGNS2(h, g, argv1.c_str(), cgnsBCMap, meshData);
  m->verify();
  //
  m->writeNative(argv2.c_str());
  // so we can see the result
  const std::string path = argv1.c_str();
  const std::size_t found = path.find_last_of("/\\");
  std::size_t index = found + 1;
  if (found == std::string::npos)
    index = 0;

  //std::cout << "FOUND " << found << " " << path.size() << " " << path << " " << index << std::endl;

#ifndef NDEBUG // debug settings, cmake double negative....
  const auto prefix = path.substr(index) + "_debug" + post;
#else // optimised setting
  const auto prefix = path.substr(index) + "_release" + post;
#endif

  const auto dim = m->getDimension();
  if (dim == 3)
  {
    {
      const auto name = prefix + "_" + std::to_string(m->getPCU()->Peers()) + "procs" + std::string("_toVTK_cellMesh");
      apf::writeVtkFiles(name.c_str(), m, 3);
    }
    {
      const auto name = prefix + "_" + std::to_string(m->getPCU()->Peers()) + "procs" + std::string("_toVTK_faceMesh");
      apf::writeVtkFiles(name.c_str(), m, 2);
    }
    {
      const auto name = prefix + "_" + std::to_string(m->getPCU()->Peers()) + "procs" + std::string("_toVTK_edgeMesh");
      apf::writeVtkFiles(name.c_str(), m, 1);
    }
    {
      const auto name = prefix + "_" + std::to_string(m->getPCU()->Peers()) + "procs" + std::string("_toVTK_vertexMesh");
      apf::writeVtkFiles(name.c_str(), m, 0);
    }
  }
  else if (dim == 2)
  {
    {
      const auto name = prefix + "_" + std::to_string(m->getPCU()->Peers()) + "procs" + std::string("_toVTK_cellMesh");
      apf::writeVtkFiles(name.c_str(), m, 2);
    }
    {
      const auto name = prefix + "_" + std::to_string(m->getPCU()->Peers()) + "procs" + std::string("_toVTK_edgeMesh");
      apf::writeVtkFiles(name.c_str(), m, 1);
    }
    {
      const auto name = prefix + "_" + std::to_string(m->getPCU()->Peers()) + "procs" + std::string("_toVTK_vertexMesh");
      apf::writeVtkFiles(name.c_str(), m, 0);
    }
  }
  else if (dim == 1)
  {
    {
      const auto name = prefix + "_" + std::to_string(m->getPCU()->Peers()) + "procs" + std::string("_toVTK_cellMesh");
      apf::writeVtkFiles(name.c_str(), m, 1);
    }
    {
      const auto name = prefix + "_" + std::to_string(m->getPCU()->Peers()) + "procs" + std::string("_toVTK_vertexMesh");
      apf::writeVtkFiles(name.c_str(), m, 0);
    }
  }

  std::string cgnsOutputName;
  if (writeCGNS)
  {
    cgnsOutputName = prefix + "_" + std::to_string(m->getPCU()->Peers()) + "procs" + "_outputFile.cgns";
    apf::writeCGNS(cgnsOutputName.c_str(), m, cgnsBCMap);
  }

  // main purpose is to call additional tests through the test harness testing.cmake
  std::function<void()> cleanUp;
  if (additionalTests)
  {
    cleanUp = additional(prefix, g, m);
    //
  }

  if (additionalTests)
  {
    // add dummy matrix to mesh
    const auto addMatrix = [](apf::Mesh2 *mesh, const int &dim) {
      apf::Field *field = nullptr;
      const std::string name = "DummyMatrix_" + std::to_string(dim);
      apf::MeshTag *tag = nullptr;
      if (dim != 0)
      {
        field = apf::createField(mesh, name.c_str(), apf::MATRIX, apf::getConstant(dim));
        if (dim == mesh->getDimension())
          tag = mesh->findTag("origCGNSGlobalElemID");
      }
      else if (dim == 0)
      {
        field = apf::createFieldOn(mesh, name.c_str(), apf::MATRIX);
        tag = mesh->findTag("origCGNSGlobalVertID");
      }
      //std::cout << "*************************** "<< dim << std::endl;
      apf::MeshIterator *it = mesh->begin(dim);
      apf::MeshEntity *elm = nullptr;
      int vals[1];
      while ((elm = mesh->iterate(it)))
      {
        apf::Matrix3x3 m(11, 12, 13, 21, 22, 23, 31, 32, 33);
        if (tag)
        {
          mesh->getIntTag(elm, tag, vals);
          m[0][0] = vals[0];
          m[1][1] = -vals[0];
          m[2][2] = vals[0] * vals[0];
        }
        apf::setMatrix(field, elm, 0, m);
      }
      mesh->end(it);
      //apf::destroyField(field);
    };

    for (int i = 0; i <= m->getDimension(); i++)
      addMatrix(m, i);

    // add dummy vector to mesh
    const auto addVector = [](apf::Mesh2 *mesh, const int &dim) {
      apf::Field *field = nullptr;
      const std::string name = "DummyVector_" + std::to_string(dim);
      apf::MeshTag *tag = nullptr;
      if (dim != 0)
      {
        field = apf::createField(mesh, name.c_str(), apf::VECTOR, apf::getConstant(dim));
        if (dim == mesh->getDimension())
          tag = mesh->findTag("origCGNSGlobalElemID");
      }
      else if (dim == 0)
      {
        field = apf::createFieldOn(mesh, name.c_str(), apf::VECTOR);
        tag = mesh->findTag("origCGNSGlobalVertID");
      }
      //std::cout << "*************************** "<< dim << std::endl;
      apf::MeshIterator *it = mesh->begin(dim);
      apf::MeshEntity *elm = nullptr;
      int vals[1];
      while ((elm = mesh->iterate(it)))
      {
        apf::Vector3 v;
        if (tag)
        {
          mesh->getIntTag(elm, tag, vals);
          v[0] = vals[0];
          v[1] = -vals[0];
          v[2] = vals[0] * vals[0];
        }
        else
        {
          v[0] = 1.0;
          v[1] = 2.0;
          v[2] = 3.0;
        }

        apf::setVector(field, elm, 0, v);
      }
      mesh->end(it);
      //apf::destroyField(field);
    };

    for (int i = 0; i <= m->getDimension(); i++)
      addVector(m, i);

    // add dummy scalar to mesh
    const auto addScalar = [](apf::Mesh2 *mesh, const int &dim) {
      apf::Field *field = nullptr;
      const std::string name = "DummyScalar_" + std::to_string(dim);
      apf::MeshTag *tag = nullptr;
      if (dim != 0)
      {
        field = apf::createField(mesh, name.c_str(), apf::SCALAR, apf::getConstant(dim));
        if (dim == mesh->getDimension())
          tag = mesh->findTag("origCGNSGlobalElemID");
      }
      else if (dim == 0)
      {
        field = apf::createFieldOn(mesh, name.c_str(), apf::SCALAR);
        tag = mesh->findTag("origCGNSGlobalVertID");
      }
      //std::cout << "*************************** "<< dim << std::endl;
      apf::MeshIterator *it = mesh->begin(dim);
      apf::MeshEntity *elm = nullptr;
      int vals[1];
      vals[0] = -1234567;
      while ((elm = mesh->iterate(it)))
      {
        double v = 1.0;
        if (tag)
        {
          mesh->getIntTag(elm, tag, vals);
          v = double(vals[0]);
        }
        apf::setScalar(field, elm, 0, v);
      }
      mesh->end(it);
      //apf::destroyField(field);
    };

    for (int i = 0; i <= m->getDimension(); i++)
      addScalar(m, i);

    if (dim == 3)
    {
      {
        const auto name = prefix + "_" + std::to_string(m->getPCU()->Peers()) + "procs" + std::string("_toVTK_cellMesh_additional");
        apf::writeVtkFiles(name.c_str(), m, 3);
      }
      {
        const auto name = prefix + "_" + std::to_string(m->getPCU()->Peers()) + "procs" + std::string("_toVTK_faceMesh_additional");
        apf::writeVtkFiles(name.c_str(), m, 2);
      }
      {
        const auto name = prefix + "_" + std::to_string(m->getPCU()->Peers()) + "procs" + std::string("_toVTK_edgeMesh_additional");
        apf::writeVtkFiles(name.c_str(), m, 1);
      }
      {
        const auto name = prefix + "_" + std::to_string(m->getPCU()->Peers()) + "procs" + std::string("_toVTK_vertexMesh_additional");
        apf::writeVtkFiles(name.c_str(), m, 0);
      }
    }
    else if (dim == 2)
    {
      {
        const auto name = prefix + "_" + std::to_string(m->getPCU()->Peers()) + "procs" + std::string("_toVTK_cellMesh_additional");
        apf::writeVtkFiles(name.c_str(), m, 2);
      }
      {
        const auto name = prefix + "_" + std::to_string(m->getPCU()->Peers()) + "procs" + std::string("_toVTK_edgeMesh_additional");
        apf::writeVtkFiles(name.c_str(), m, 1);
      }
      {
        const auto name = prefix + "_" + std::to_string(m->getPCU()->Peers()) + "procs" + std::string("_toVTK_vertexMesh_additional");
        apf::writeVtkFiles(name.c_str(), m, 0);
      }
    }
    else if (dim == 1)
    {
      {
        const auto name = prefix + "_" + std::to_string(m->getPCU()->Peers()) + "procs" + std::string("_toVTK_cellMesh_additional");
        apf::writeVtkFiles(name.c_str(), m, 1);
      }
      {
        const auto name = prefix + "_" + std::to_string(m->getPCU()->Peers()) + "procs" + std::string("_toVTK_vertexMesh_additional");
        apf::writeVtkFiles(name.c_str(), m, 0);
      }
    }

    if (writeCGNS)
    {
      // what this one to be re-read if doing re-reading so that the dummy variables (vector/matrix/scalar) are there
      cgnsOutputName = prefix + "_" + std::to_string(m->getPCU()->Peers()) + "procs" + "_additional_outputFile.cgns";
      apf::writeCGNS(cgnsOutputName.c_str(), m, cgnsBCMap);
    }
    //
    {
      apf::Field *field = nullptr;
      const std::string name = "displace";
      field = apf::createFieldOn(m, name.c_str(), apf::VECTOR);
      apf::MeshIterator *it = m->begin(0);
      apf::MeshEntity *elm = nullptr;
      while ((elm = m->iterate(it)))
      {
        apf::Vector3 v;
        v[0] = 0.25;
        v[1] = 0.25;
        v[2] = 0.35;
        apf::setVector(field, elm, 0, v);
      }
      m->end(it);
      apf::displaceMesh(m, field);
      if (writeCGNS)
      {
        std::string cgnsOutputName = prefix + "_" + std::to_string(m->getPCU()->Peers()) + "procs" + "_additional_moved_outputFile.cgns";
        apf::writeCGNS(cgnsOutputName.c_str(), m, cgnsBCMap);
      }
    }
  }

  if (additionalTests)
  {
    //cleanUp();
  }
  else if (!additionalTests)
  {
    m->destroyNative();
    apf::destroyMesh(m);
  }
  //
  return cgnsOutputName;
}

int main(int argc, char **argv)
{
#ifdef HAVE_CGNS
  MPI_Init(&argc, &argv);
  {
  pcu::PCU PCUObj = pcu::PCU(MPI_COMM_WORLD);
  pumi_load_pcu(&PCUObj);
  lion_set_verbosity(1);
  bool additionalTests = false;
  if (argc < 3)
  {
    if (!PCUObj.Self())
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
      if (!PCUObj.Self())
        printf("Usage: %s <in .cgns> <out .smb> additional\n", argv[0]);
      MPI_Finalize();
      exit(EXIT_FAILURE);
      return -1;
    }
  }
  else if (argc > 4)
  {
    if (!PCUObj.Self())
      printf("Usage: %s <in .cgns> <out .smb>\n", argv[0]);
    MPI_Finalize();
    exit(EXIT_FAILURE);
    return -1;
  }

  // Phase 1
  std::string cgnsOutputName;
  {
    apf::CGNSBCMap cgnsBCMap;
    std::vector<std::pair<std::string, std::string>> meshData;
    cgnsOutputName = doit(cgnsBCMap, argv[1], argv[2], "", additionalTests, additionalTests, meshData, &PCUObj);
  }
  // Phase 2
  if (additionalTests)
  {
    std::vector<std::pair<std::string, std::string>> meshData;
    meshData.push_back(std::make_pair("Vertex", "DummyScalar_0"));
    meshData.push_back(std::make_pair("Vertex", "DummyVector_0"));
    meshData.push_back(std::make_pair("Vertex", "DummyMatrix_0"));

    meshData.push_back(std::make_pair("CellCenter", "DummyScalar_3"));
    meshData.push_back(std::make_pair("CellCenter", "DummyVector_3"));
    meshData.push_back(std::make_pair("CellCenter", "DummyMatrix_3"));

    meshData.push_back(std::make_pair("CellCenter", "DummyScalar_2"));
    meshData.push_back(std::make_pair("CellCenter", "DummyVector_2"));
    meshData.push_back(std::make_pair("CellCenter", "DummyMatrix_2"));

    meshData.push_back(std::make_pair("CellCenter", "DummyScalar_1"));
    meshData.push_back(std::make_pair("CellCenter", "DummyVector_1"));
    meshData.push_back(std::make_pair("CellCenter", "DummyMatrix_1"));
    apf::CGNSBCMap cgnsBCMap;
    std::cout << "RE-READING " << std::endl;
    doit(cgnsBCMap, cgnsOutputName.c_str(), "tempy.smb", "_reread", false, true, meshData, &PCUObj);
  }
  //
  }
  MPI_Finalize();
  return 0;
#else
  PCU_ALWAYS_ASSERT_VERBOSE(true == false,
                            "Build with ENABLE_CGNS to allow this functionality.");
  exit(EXIT_FAILURE);
  return -1;
#endif
}
