#include <lionPrint.h>
#include <gmi_null.h>
#include <apfMDS.h>
#include <apf.h>
#include <apfMesh2.h>
#include <apfNumbering.h>
#include <pcu_util.h>
#include <cstdlib>

void createMesh(gmi_model*& g, apf::Mesh2*& m, int n, pcu::PCU *PCUObj)
{
  gmi_register_null();
  g = gmi_load(".null");
  m = apf::makeEmptyMdsMesh(g, 1, false, PCUObj);
  apf::ModelEntity* left = m->findModelEntity(0, 0);
  apf::ModelEntity* right = m->findModelEntity(0, 1);
  apf::ModelEntity* mid = m->findModelEntity(1, 1);
  apf::NewArray<apf::MeshEntity*> v(n);
  for (int i = 0; i < n; ++i) {
    apf::ModelEntity* c;
    if (i == 0)
      c = left;
    else if (i == n-1)
      c = right;
    else
      c = mid;
    v[i] = m->createVert(c);
    double x = double(i)/double(n-1);
    m->setPoint(v[i],0,apf::Vector3(x,0,0));
  }
  for (int i = 0; i < n-1; ++i) {
    apf::MeshEntity* ev[2];
    ev[0] = v[i];
    ev[1] = v[i+1];
    m->createEntity(apf::Mesh::EDGE, mid, ev);
  }
  m->acceptChanges();
  m->verify();
}

void test(apf::Mesh2* m)
{
  /*
  apf::GlobalNumbering* gn = 
    apf::createGlobalNumbering(
        m,"nodes",apf::getShape(m->getCoordinateField()));
   */
  int dim = 0;
  int tag = 0;
  apf::ModelEntity* v = m->findModelEntity(dim,tag);
  apf::DynamicArray<apf::Node> nodes;
  apf::getNodesOnClosure(m,v,nodes);
  printf("%lu nodes on geometric entity with dim: "
         "%d and tag: %d\n",nodes.getSize(),dim,tag);
  /*
  for (std::size_t i=0; i < nodes.getSize(); ++i)
    printf("node has id: %lu\n",i);
   */
}

int main(int argc, char** argv)
{
  /** 1 - number of points
    * 2 - model output name
    * 3 = mesh output name **/
  PCU_ALWAYS_ASSERT(argc==4);
  pcu::Init(&argc,&argv);
  {
  pcu::PCU pcu_obj;
  lion_set_verbosity(1);
  gmi_model* g;
  apf::Mesh2* m;
  int nverts = atoi(argv[1]);
  PCU_ALWAYS_ASSERT(2 <= nverts);
  PCU_ALWAYS_ASSERT(nverts <= 1000);
  createMesh(g,m,nverts,&pcu_obj);
  test(m);
  gmi_write_dmg(g,argv[2]);
  m->writeNative(argv[3]);
  m->destroyNative();
  apf::destroyMesh(m);
  }
  pcu::Finalize();
}
