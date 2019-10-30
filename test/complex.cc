#include <apf.h>
#include <apfComplex.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <apfShape.h>
#include <gmi_mesh.h>
#include <PCU.h>
#include <lionPrint.h>
#include <pcu_util.h>

int main(int ac, char** av)
{
  // requires C++ complex, not C complex
  using namespace std::complex_literals;
  PCU_ALWAYS_ASSERT(ac == 3);
  MPI_Init(&ac,&av);
  PCU_Comm_Init();
  gmi_register_mesh();
  apf::Mesh2* m = apf::loadMdsMesh(av[1], av[2]);
  static const int cmps = 2;
  apf::ComplexField* f = apf::createComplexPackedField(m,"foo",cmps);
  int vrts = 0;
  apf::zeroField(f);
  double_complex zero = {0,0};
  double_complex r[cmps];
  double_complex w[cmps] = { {1,2},{2,1} };
  {
    apf::MeshEntity* vert = NULL;
    apf::MeshIterator* it = m->begin(0);
    while ((vert = m->iterate(it)))
    {
      apf::getComponents(f,vert,0,&r[0]);
      PCU_ALWAYS_ASSERT(r[0] == zero && r[1] == zero);
      apf::setComponents(f,vert,0,&w[0]);
      apf::getComponents(f,vert,0,&r[0]);
      PCU_ALWAYS_ASSERT(r[0] == w[0] && r[1] == w[1]);
      ++vrts;
    }
    m->end(it);
  }
  apf::freeze(f);
  PCU_ALWAYS_ASSERT(apf::isFrozen(f));
  {
    apf::MeshEntity* vert = NULL;
    apf::MeshIterator* it = m->begin(0);
    while ((vert = m->iterate(it)))
    {
      apf::getComponents(f,vert,0,&r[0]);
      PCU_ALWAYS_ASSERT(r[0] == w[0] && r[1] == w[1]);
    }
    m->end(it);
  }
  double_complex * foo_arr = apf::getComplexArrayData(f);
  for(int vv = 0; vv < vrts; ++vv)
    PCU_ALWAYS_ASSERT(foo_arr[cmps*vv] == r[0] && foo_arr[cmps*vv+1] == r[1]);
  apf::unfreeze(f);
  foo_arr = NULL; // im certain if we tried to access this we would get an error,
                  // but there isn't a good way to prevent this without changing the API
  PCU_ALWAYS_ASSERT(!apf::isFrozen(f));

  // TODO: test getComponents at a parametric location
  //       getshapevalues and getshapegrads only really use the meshelement so dont bother
  {
    apf::MeshEntity* e = NULL;
    int dim = m->getDimension();
    apf::MeshIterator* it = m->begin(dim);
    while ((e = m->iterate(it)))
    {
      apf::MeshElement * me = apf::createMeshElement(m,e);
      apf::ComplexElement * el1 = apf::createElement(f,me);
      apf::ComplexElement * el2 = apf::createElement(f,e);
      apf::MeshElement * me1 = apf::getMeshElement(el1);
      apf::MeshElement * me2 = apf::getMeshElement(el2);
      PCU_ALWAYS_ASSERT(me1 == me && me2 == NULL);
      apf::MeshEntity * e1 = apf::getMeshEntity(el1);
      apf::MeshEntity * e2 = apf::getMeshEntity(el2);
      PCU_ALWAYS_ASSERT(e1 == e && e2 == e);
      int nds1 = apf::countNodes(el1);
      int nds2 = apf::countNodes(el2);
      PCU_ALWAYS_ASSERT(nds1 == nds2);
      apf::NewArray<double_complex> dofs1(cmps*nds1);
      apf::NewArray<double_complex> dofs2(cmps*nds2);
      apf::getPackedNodes(el1,dofs1);
      apf::getPackedNodes(el2,dofs2);
      for(int nd = 0; nd < nds1; ++nd)
      {
        PCU_ALWAYS_ASSERT(dofs1[nd*cmps] == w[0] && dofs1[nd*cmps+1] == w[1]);
        PCU_ALWAYS_ASSERT(dofs2[nd*cmps] == w[0] && dofs2[nd*cmps+1] == w[1]);
      }
      apf::destroyElement(el2);
      apf::destroyElement(el1);
      apf::destroyMeshElement(me);
    }
    m->end(it);
  }
  m->destroyNative();
  apf::destroyMesh(m);
  PCU_Comm_Free();
  MPI_Finalize();
}
