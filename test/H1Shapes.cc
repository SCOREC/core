#include <PCU.h>
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


// User defined vector functions E(x,y,z) of order up to 6
void E_exact(const apf::Vector3& x, apf::Vector3& value, int p);

void testH1(
    apf::Mesh2* m,
    const apf::Vector3& testXi,
    int ndOrder, int exactOrder);

int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();

  lion_set_verbosity(1);

  if (argc != 3) {
    if(0==PCU_Comm_Self())
      std::cerr << "usage: " << argv[0]
        << " <model.dmg or .null> <mesh.smb>\n";
    return EXIT_FAILURE;
  }

  gmi_register_mesh();
  gmi_register_null();

  gmi_model* g = gmi_load(argv[1]);
  apf::Mesh2* m = apf::loadMdsMesh(g,argv[2]);
  m->verify();

  for (int i = 1; i <= 6; i++) {
    testH1(
    	m, /* mesh */
	apf::Vector3(1./7., 1./11., 1./3.), /* test point */
	i, /* order of h1 field */
	i); /* order of test field */
  }

  apf::destroyMesh(m);
  PCU_Comm_Free();
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

void testH1(
    apf::Mesh2* m,
    const apf::Vector3& testXi,
    int ndOrder, int exactOrder)
{
  // Loop over all nodes and set scalar dofs.
  int dim = m->getDimension();
  apf::Field* h1Field = apf::createField(
      m, "h1_test", apf::VECTOR, apf::getH1Shape(ndOrder));
  apf::MeshEntity* ent;
  apf::MeshIterator* it;

  for (int d = 0; d <= dim; d++) {
    if (!h1Field->getShape()->countNodesOn(apf::Mesh::simplexTypes[d])) {
      if(0==PCU_Comm_Self())
        lion_oprint(1, "no nodes in dimension %d\n", d);
      continue;
    }
    else
      if(0==PCU_Comm_Self())
        lion_oprint(1, "computing dofs for dimension %d\n", d);
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


  // Verify that interpolated solution field agrees with exact field.
  for (int d = 0; d <= dim; d++) {

    double L2ErrorE = 0.;
    it = m->begin(d);
    int count = 0;
    while( (ent = m->iterate(it)) ) {
      apf::MeshElement* me = apf::createMeshElement(m, ent);
      apf::Vector3 x;
      apf::mapLocalToGlobal(me, testXi, x);
      apf::Vector3 eFieldExact;
      E_exact(x, eFieldExact, exactOrder);

      // obtain interpolated value
      apf::Element* el = apf::createElement(h1Field, me);
      apf::Vector3 eFieldValue;
      apf::getVector(el, testXi, eFieldValue);

      double err = ((eFieldValue - eFieldExact) * (eFieldValue - eFieldExact));
      err /= (eFieldExact * eFieldExact); // normalization factor
      L2ErrorE += err;
      apf::destroyMeshElement(me);
      apf::destroyElement(el);
      count++;
    }
    m->end(it);

    // check for field interpolation
    if(0==PCU_Comm_Self())
      lion_oprint(1, "L2ErrorE for entities of dimension %d is %e\n", d, L2ErrorE);
    PCU_ALWAYS_ASSERT_VERBOSE(L2ErrorE < 1.e-12,
	"Fields were not interpolated correctly!");

  }
  m->removeField(h1Field);
  apf::destroyField(h1Field);
}
