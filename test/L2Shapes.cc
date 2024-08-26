#include <gmi_mesh.h>
#include <gmi_null.h>
#include <apf.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <apfShape.h>
#include <apfField.h>
#include <lionPrint.h>
#include <pcu_util.h>

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <string>
#include <math.h>

using namespace std;

// User defined vector functions E(x,y,z) of order up to 6
void E_exact(const apf::Vector3& x, apf::Vector3& value, int p);

void testL2(
    apf::Mesh2* m,
    const apf::Vector3& testXi,
    int ndOrder, int exactOrder);

void testL2writeNative(
    apf::Mesh2* m,
    int ndOrder, int exactOrder);

int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  {
  pcu::PCU PCUObj = pcu::PCU(MPI_COMM_WORLD);

  lion_set_verbosity(1);

  if (argc != 3) {
    if(0==PCUObj.Self())
      std::cerr << "usage: " << argv[0]
        << " <model.dmg or .null> <mesh.smb>\n";
    return EXIT_FAILURE;
  }

  gmi_register_mesh();
  gmi_register_null();

  gmi_model* g = gmi_load(argv[1]);
  apf::Mesh2* m = apf::loadMdsMesh(g,argv[2],&PCUObj);
  m->verify();

  for (int i = 0; i <= 6; i++) {
    testL2(
    	m, /* mesh */
	apf::Vector3(1./7., 1./11., 1./3.), /* test point */
	i, /* order of l2 field */
	i); /* order of test field */
  }

  testL2writeNative(m, 3, 3);

  apf::destroyMesh(m);
  }
  MPI_Finalize();
}

void E_exact(const apf::Vector3& x, apf::Vector3& value, int p)
{
  // Polynomial coefficients for each component of exact vector field
  double a[7] = { 1.0, -1.0,  2., -2., -1.0,  1.0, -1.0};
  double b[7] = {-2.0,  1.0, -2.,  2., -1.0, -1.0,  1.0};
  double c[7] = { 3.0,  0.0, -1.,  1., -1.0,  1.0,  0.0};

  value[0] = 0.0;
  value[1] = 0.0;
  value[2] = 0.0;
  for (int i = p;  i >= 0; i--) {
    value[0] += pow(x[0],p)*a[p];
    value[1] += pow(x[1],p)*b[p];
    value[2] += pow(x[2],p)*c[p];
  }
}

void testL2(
    apf::Mesh2* m,
    const apf::Vector3& testXi,
    int ndOrder, int exactOrder)
{
  // Loop over all nodes and set scalar dofs.
  int dim = m->getDimension();
  apf::Field* l2Field = apf::createField(
      m, "l2_test", apf::VECTOR, apf::getL2Shape(ndOrder, apf::Mesh::simplexTypes[dim]));
  apf::MeshEntity* ent;
  apf::MeshIterator* it;

  for (int d = 0; d <= dim; d++) {
    if (!l2Field->getShape()->countNodesOn(apf::Mesh::simplexTypes[d])) {
      if(0==m->getPCU()->Self())
        lion_oprint(1, "no nodes in dimension %d\n", d);
      continue;
    }
    else
      if(0==m->getPCU()->Self())
        lion_oprint(1, "computing dofs for dimension %d\n", d);
    it = m->begin(d);
    while( (ent = m->iterate(it)) ) {
      int type = m->getType(ent);
      int non = l2Field->getShape()->countNodesOn(type);
      apf::MeshElement* me = apf::createMeshElement(m, ent);
      for (int i = 0; i < non; i++)
      {
        apf::Vector3 xi, p, value;
        l2Field->getShape()->getNodeXi(type, i, xi);
        apf::mapLocalToGlobal(me, xi, p);
        E_exact(p, value, exactOrder);

      	apf::setVector(l2Field, ent, i, value);
      }
      apf::destroyMeshElement(me);
    }
    m->end(it);
  }


  // Verify that interpolated solution field agrees with exact field.
  double L2ErrorE = 0.;
  it = m->begin(dim);
  while( (ent = m->iterate(it)) ) {
    apf::MeshElement* me = apf::createMeshElement(m, ent);
    apf::Vector3 x;
    apf::mapLocalToGlobal(me, testXi, x);
    apf::Vector3 eFieldExact;
    E_exact(x, eFieldExact, exactOrder);

    // obtain interpolated value
    apf::Element* el = apf::createElement(l2Field, me);
    apf::Vector3 eFieldValue;
    apf::getVector(el, testXi, eFieldValue);

    double err = ((eFieldValue - eFieldExact) * (eFieldValue - eFieldExact));
    err /= (eFieldExact * eFieldExact); // normalization factor
    L2ErrorE += err;
    apf::destroyMeshElement(me);
    apf::destroyElement(el);
  }
  m->end(it);

  // check for field interpolation
  if(0==m->getPCU()->Self())
    lion_oprint(1, "L2ErrorE is %e\n", L2ErrorE);
  PCU_ALWAYS_ASSERT_VERBOSE(L2ErrorE < 1.e-16,
      "Fields were not interpolated correctly!");

  m->removeField(l2Field);
  apf::destroyField(l2Field);
}

void testL2writeNative(
    apf::Mesh2* m,
    int ndOrder, int exactOrder)
{
  // Loop over all nodes and set scalar dofs.
  int dim = m->getDimension();
  apf::Field* l2Field = apf::createField(
      m, "l2_test", apf::VECTOR, apf::getL2Shape(ndOrder, apf::Mesh::simplexTypes[dim]));
  apf::MeshEntity* ent;
  apf::MeshIterator* it;

  for (int d = 0; d <= dim; d++) {
    if (!l2Field->getShape()->countNodesOn(apf::Mesh::simplexTypes[d])) {
      if(0==m->getPCU()->Self())
        lion_oprint(1, "no nodes in dimension %d\n", d);
      continue;
    }
    else
      if(0==m->getPCU()->Self())
        lion_oprint(1, "computing dofs for dimension %d\n", d);
    it = m->begin(d);
    while( (ent = m->iterate(it)) ) {
      int type = m->getType(ent);
      int non = l2Field->getShape()->countNodesOn(type);
      apf::MeshElement* me = apf::createMeshElement(m, ent);
      for (int i = 0; i < non; i++)
      {
        apf::Vector3 xi, p, value;
        l2Field->getShape()->getNodeXi(type, i, xi);
        apf::mapLocalToGlobal(me, xi, p);
        E_exact(p, value, exactOrder);
        apf::setVector(l2Field, ent, i, value);
      }
      apf::destroyMeshElement(me);
    }
    m->end(it);
  }

  // project to nodal field
  apf::Field* projectedL2Field = apf::createField(m, "projected_l2_test", apf::VECTOR,
      apf::getLagrange(1));
  apf::projectL2Field(projectedL2Field, l2Field);

  // write to native
  // 1- write the mesh to native
  // 2- read the mesh back in and make sure fields are on the mesh
  // 3- clean up the newly loaded mesh
  m->writeNative("L2Shape_test_mesh.smb");
  apf::Mesh2* m2 = apf::loadMdsMesh(".null", "./L2Shape_test_mesh.smb", m->getPCU());
  int fCount = 0;
  for (int i = 0; i < m2->countFields(); i++) {
    if(0==m2->getPCU()->Self())
      lion_oprint(1, "field %d's name and shape are %s and %s\n", i,
	  m2->getField(i)->getName(), m2->getField(i)->getShape()->getName());
    fCount++;
  }
  PCU_ALWAYS_ASSERT_VERBOSE(fCount == 2, "Expecting 2 fields on m2 at this point");

  apf::destroyMesh(m2);
  // clean up the fields
  apf::destroyField(l2Field);
  apf::destroyField(projectedL2Field);
}
