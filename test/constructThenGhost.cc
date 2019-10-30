#include <gmi_mesh.h>
#include <gmi_null.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <apfConvert.h>
#include <apf.h>
#include <PCU.h>
#include <lionPrint.h>
#include <pcu_util.h>
#include <pumi.h>
#include <algorithm>

int main(int argc, char** argv)
{
  PCU_ALWAYS_ASSERT(argc==3);
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  lion_set_verbosity(1);
  gmi_register_mesh();
  gmi_register_null();
  int* conn;
  double* coords;
  int nelem;
  int etype;
  int nverts;

  apf::Mesh2* m = apf::loadMdsMesh(argv[1],argv[2]);
  int dim = m->getDimension();
  extractCoords(m, coords, nverts);
  destruct(m, conn, nelem, etype);
  m->destroyNative();
  apf::destroyMesh(m);

  gmi_model* model = gmi_load(".null");
  m = apf::makeEmptyMdsMesh(model, dim, false);
  apf::GlobalToVert outMap;
  apf::construct(m, conn, nelem, etype, outMap);
  delete [] conn;
  apf::alignMdsRemotes(m);
  apf::deriveMdsModel(m);
  apf::setCoords(m, coords, nverts, outMap);
  delete [] coords;
  outMap.clear();
  m->verify();

  if (!pumi_rank()) std::cout<<"model/mesh converted to pumi instance\n";

  //create the pumi instance to use pumi api's
  pGeom g=pumi_geom_load(model);
  pMesh pm = pumi_mesh_load(m);
  pumi_mesh_verify(pm);

  //create an element field
  const int mdim = pumi_mesh_getDim(pm);
  pShape s = pumi_shape_getConstant(mdim);
  const int dofPerElm = 1;
  pField f = pumi_field_create(pm, "elmField", dofPerElm, PUMI_PACKED, s);

  pMeshIter it = pm->begin(mdim);
  pMeshEnt e;
  double v = 0;
  while ((e = m->iterate(it)))
    pumi_node_setField(f,e,0,&v);
  pm->end(it);

  const int ghost = mdim;
  const int bridge = ghost-1;
  const int numLayers = 1;
  const int ghostAcrossCopies = 1;
  pumi_ghost_createLayer(pm,bridge,ghost,numLayers,ghostAcrossCopies);

  it = pm->begin(mdim);
  v = 1;
  while ((e = m->iterate(it))) {
    if (!pumi_ment_isGhost(e))
      pumi_node_setField(f,e,0,&v);
  }
  pm->end(it);

  // only the owned elements will have a elmField value of 1
  pumi_mesh_write(pm, "beforeSync", "vtk");

  pumi_field_synchronize(f);

  // owned and ghosted elements will have a elmField value of 1
  pumi_mesh_write(pm, "afterSync", "vtk");

  // clean-up 
  pumi_field_delete(f);
  pumi_ghost_delete(pm);

  pumi_geom_delete(g);
  pumi_mesh_delete(pm);
  PCU_Comm_Free();
  MPI_Finalize();
}
