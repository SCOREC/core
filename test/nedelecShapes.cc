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

void testNedelec(
    apf::Mesh2* m,
    const apf::Vector3& testXi,
    int ndOrder, int exactOrder);

int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  {
  pcu::PCU PCUObj = pcu::PCU(MPI_COMM_WORLD);

  lion_set_verbosity(0);

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


  // constant fields have to be interpolated exactly by any order Nedelec shape
  for (int order = 1; order <= 6; order++) {
    testNedelec(
	m, /* mesh */
	apf::Vector3(1./4., 1./5., 1./6.), /* test point */
	order, /* order of Nedelec field */
	0); /* constant field */
  }

  // Fields of order p are  interpolated exactly by Nedelec shapes of order p+1
  for (int i = 1; i <= 6; i++) {
    testNedelec(
	m, /* mesh */
	apf::Vector3(1./4., 1./5., 1./6.), /* test point */
	i, /* order of Nedelec field */
	i-1); /* order of test field */
  }

  apf::destroyMesh(m);
  }
  MPI_Finalize();
}

void E_exact(const apf::Vector3& x, apf::Vector3& value, int p)
{
  // Polynomial coefficients for each component of exact vector field
  double a[6] = { 1.0, -1.0,  2., -2., -1.0,  1.0};
  double b[6] = {-2.0,  1.0, -2.,  2., -1.0, -1.0};
  double c[6] = { 3.0,  0.0, -1.,  0., -1.0,  1.0};

  value[0] = 0.0;
  value[1] = 0.0;
  value[2] = 0.0;
  for (int i = p;  i >= 0; i--) {
    value[0] += pow(x[0],p)*a[p];
    value[1] += pow(x[1],p)*b[p];
    value[2] += pow(x[2],p)*c[p];
  }
}



void testNedelec(
    apf::Mesh2* m,
    const apf::Vector3& testXi,
    int ndOrder, int exactOrder)
{
  apf::Field* ndField = apf::createField(
      m, "nedelec_test", apf::SCALAR, apf::getNedelec(ndOrder));

  // Loop over all nodes and set scalar dofs.
  int dim = m->getDimension();
  apf::MeshEntity* ent;
  apf::MeshIterator* it;

  for (int d = 0; d <= dim; d++) {
    if (!ndField->getShape()->countNodesOn(apf::Mesh::simplexTypes[d])) {
      lion_oprint(1, "no nodes in dimension %d\n", d);
      continue;
    }
    else
      lion_oprint(1, "computing dofs for dimension %d\n", d);
    it = m->begin(d);
    while( (ent = m->iterate(it)) ) {
      int type = m->getType(ent);
      int non = ndField->getShape()->countNodesOn(type);
      apf::MeshElement* me = apf::createMeshElement(m, ent);
      for (int i = 0; i < non; i++)
      {
        apf::Vector3 xi, p, value;
        ndField->getShape()->getNodeXi(type, i, xi);
        apf::mapLocalToGlobal(me, xi, p);
        E_exact(p, value, exactOrder);

        apf::Matrix3x3 J;
        apf::getJacobian( me, xi, J);

        apf::Vector3 t;
        ndField->getShape()->getNodeTangent(type, i, t);


	// dof is t . J^T value
	// note getJacobian returns the transpose of J, so no need to transpose
	// J again.
        apf::Vector3 temp = J * value;
        double dof = temp * t;
      	apf::setScalar(ndField, ent, i, dof);
      }
      apf::destroyMeshElement(me);
    }
    m->end(it);
  }


  // Verify that interpolated solution field agrees with exact field.
  double L2ErrorE = 0.;
  double L2ErrorCurlE = 0.;
  it = m->begin(3);
  while( (ent = m->iterate(it)) ) {
    apf::MeshElement* me = apf::createMeshElement(m, ent);
    apf::Vector3 x;
    apf::mapLocalToGlobal(me, testXi, x);
    apf::Vector3 eFieldExact;
    E_exact(x, eFieldExact, exactOrder);

    // obtain interpolated value
    apf::Element* el = apf::createElement(ndField, me);
    apf::Vector3 eFieldValue;
    apf::getVector(el, testXi, eFieldValue);
    // obtain curl field values
    apf::Vector3 eCurlValue;
    apf::getCurl(el, testXi, eCurlValue);

    L2ErrorE += ((eFieldValue - eFieldExact) * (eFieldValue - eFieldExact))
      / (eFieldExact * eFieldExact); // normalization factor
    L2ErrorCurlE += eCurlValue * eCurlValue;
    apf::destroyMeshElement(me);
    apf::destroyElement(el);
  }
  m->end(it);

  // check for field interpolation
  PCU_ALWAYS_ASSERT_VERBOSE(L2ErrorE < 1.e-16,
      "Fields were not interpolated correctly!");
  // check for curl (only applies for constant field)
  if (exactOrder == 0)
    PCU_ALWAYS_ASSERT_VERBOSE(L2ErrorCurlE < 1.e-16,
	"Curls are expected to be zero!");

  apf::destroyField(ndField);
}
