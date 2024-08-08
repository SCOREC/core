#include <apf.h>
#include <apfMDS.h>
#include <apfShape.h>
#include <apfField.h>
#include <gmi_mesh.h>
#include <ma.h>
#include <maShape.h>
#include <crv.h>
#include <lionPrint.h>
#ifdef HAVE_SIMMETRIX
#include <SimUtil.h>
#include <MeshSim.h>
#include <gmi_sim.h>
#include <SimModel.h>
#endif
#include <cassert>
#include <stdlib.h>
#include <iostream>
#include <memory>

void E_exact(const apf::Vector3& x, apf::Vector3& value, int p);

apf::Field* addH1Field(
    apf::Mesh2* m,
    int ndOrder, int exactOrder);

double testH1Field(
    apf::Mesh2* m,
    apf::Field* f,
    const apf::Vector3& testXi,
    int exactOrder);

void testCurveAdapt(
    const char* modelFile,
    const char* meshFile,
    pcu::PCU *PCUObj,
    const int mesh_order,
    const int exact_order,
    const int field_order);

int main(int argc, char** argv)
{
  const char* modelFile = argv[1];
  const char* meshFile = argv[2];
  MPI_Init(&argc,&argv);
  {
  auto PCUObj = std::unique_ptr<pcu::PCU>(new pcu::PCU(MPI_COMM_WORLD));
  lion_set_verbosity(1);
#ifdef HAVE_SIMMETRIX
  MS_init();
  SimModel_start();
  Sim_readLicenseFile(0);
  gmi_sim_start();
  gmi_register_sim();
#endif
  gmi_register_mesh();


  /* Note on choices of the orders:
   * -- mesh geometry if of order p   [that is x(xi) is a polynomial of order p in xi]
   * -- test (exact) field is order q [that is F(x) is a polynomial of order q in x]
   * -- then the finite element basis has to be of order at least p*q. This is because
   *    F(x) = F(x(xi)) will be a polynomial of order p*q in xi.
   */

  // linear adapt
  testCurveAdapt(modelFile, meshFile, PCUObj.get(),
      1 /*mesh_order*/,
      2 /*exact_order*/,
      2 /*field_order*/);

  // quadratic adapts
  testCurveAdapt(modelFile, meshFile, PCUObj.get(),
      2 /*mesh_order*/,
      2 /*exact_order*/,
      4 /*field_order*/);
  testCurveAdapt(modelFile, meshFile, PCUObj.get(),
      2 /*mesh_order*/,
      3 /*exact_order*/,
      6 /*field_order*/);

  // cubic adapt
  testCurveAdapt(modelFile, meshFile, PCUObj.get(),
      3 /*mesh_order*/,
      2 /*exact_order*/,
      6 /*field_order*/);

  }
#ifdef HAVE_SIMMETRIX
  gmi_sim_stop();
  Sim_unregisterAllKeys();
  MS_exit();
  SimModel_stop();
#endif
  MPI_Finalize();
}


void E_exact(const apf::Vector3& x, apf::Vector3& value, int p)
{
  // Polynomial coefficients for each component of exact vector field
  double a[7] = { 1.0, -1.0,  2., -2., -1.0,  1.0, -1.0};
  double b[7] = {-2.0,  1.0, -2.,  2., -1.0, -1.0,  1.0};
  double c[7] = { 3.0,  0.0, -1.,  100., -1.0,  1.0,  0.0};

  value[0] = 0.0;
  value[1] = 0.0;
  value[2] = 0.0;
  for (int i = p;  i >= 0; i--) {
    value[0] += pow(x[0],p)*a[p];
    value[1] += pow(x[1],p)*b[p];
    value[2] += pow(x[2],p)*c[p];
  }
}

apf::Field* addH1Field(
    apf::Mesh2* m,
    int ndOrder, int exactOrder)
{
  // Loop over all nodes and set scalar dofs.
  int dim = m->getDimension();
  apf::Field* h1Field = apf::createField(
      m, "h1_test", apf::VECTOR, apf::getH1Shape(ndOrder));
  apf::MeshEntity* ent;
  apf::MeshIterator* it;

  for (int d = 0; d <= dim; d++) {
    if (!h1Field->getShape()->countNodesOn(apf::Mesh::simplexTypes[d]))
      continue;

    it = m->begin(d);
    int count = 0;
    while( (ent = m->iterate(it)) ) {
      int type = m->getType(ent);
      int non = h1Field->getShape()->countNodesOn(type);
      apf::MeshElement* me = apf::createMeshElement(m, ent);
      for (int i = 0; i < non; i++)
      {
        apf::Vector3 xi, p, value;
        h1Field->getShape()->getNodeXi(type, i, xi);
        apf::mapLocalToGlobal(me, xi, p);
        E_exact(p, value, exactOrder);

        apf::setVector(h1Field, ent, i, value);
      }
      apf::destroyMeshElement(me);
      count++;
    }
    m->end(it);
  }

  return h1Field;
}

double testH1Field(
    apf::Mesh2* m,
    apf::Field* f,
    const apf::Vector3& testXi,
    int exactOrder)
{

  int dim = m->getDimension();
  int count;
  apf::MeshIterator* it;
  apf::MeshEntity* ent;

  // Verify that interpolated solution field agrees with exact field.
  double L2ErrorE = 0.;
  it = m->begin(dim);
  count = 0;
  while( (ent = m->iterate(it)) ) {
    apf::MeshElement* me = apf::createMeshElement(m, ent);
    apf::Vector3 x;
    apf::mapLocalToGlobal(me, testXi, x);
    apf::Vector3 eFieldExact;
    E_exact(x, eFieldExact, exactOrder);

    // obtain interpolated value
    apf::Element* el = apf::createElement(f, me);
    apf::Vector3 eFieldValue;
    apf::getVector(el, testXi, eFieldValue);

    double err = ((eFieldValue - eFieldExact) * (eFieldValue - eFieldExact));
    /* err /= (eFieldExact * eFieldExact); // normalization factor */
    L2ErrorE += err;
    apf::destroyMeshElement(me);
    apf::destroyElement(el);
    count++;
  }
  m->end(it);

  return L2ErrorE;
}

void testCurveAdapt(
    const char* modelFile,
    const char* meshFile,
    pcu::PCU *PCUObj,
    const int mesh_order,
    const int exact_order,
    const int field_order)
{

  apf::Mesh2* m = apf::loadMdsMesh(modelFile,meshFile,PCUObj);
  m->verify();

  if (mesh_order > 1) {
    crv::BezierCurver bc(m, mesh_order, 0);
    bc.run();
  }

  apf::Field* f = addH1Field(m, field_order, exact_order);
  double l2ErrorBefore = testH1Field(m, f, apf::Vector3(1./3., 1./4., 1./5.), exact_order);

  ma::Input* in = ma::makeAdvanced(ma::configureUniformRefine(m,1));
  // Snap is off for solutions transfer testing.
  in->shouldSnap = false;
  in->goodQuality = 0.3*0.3*0.3;
  in->shouldFixShape = false;
  if (mesh_order > 1)
    crv::adapt(in);
  else
    ma::adapt(in);

  double l2ErrorAfter = testH1Field(m, f, apf::Vector3(1./3., 1./4., 1./5.), exact_order);

  lion_oprint(1, "mesh/exact/fields orders  %d/%d/%d,  before/after errors %e/%e\n",
      mesh_order, exact_order, field_order, l2ErrorBefore, l2ErrorAfter);

  PCU_ALWAYS_ASSERT(l2ErrorAfter < 1.e-16);

  m->destroyNative();
  apf::destroyMesh(m);
}
