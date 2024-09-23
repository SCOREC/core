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
void M_exact(const apf::Vector3& x, apf::Vector3& value, int p);

void testH1(
    apf::Mesh2* m,
    const apf::Vector3& testXi,
    int ndOrder, int exactOrder,
    int dataType); /* apf::VECTOR or apf::MATRIX */

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

  // test fields interpolating a user-defined vector field
  for (int i = 1; i <= 6; i++) {
    if(0==PCUObj.Self())
      lion_oprint(1, "----TESTING VECTOR FIELD OF ORDER %d----\n", i);
    testH1(
        m, /* mesh */
        apf::Vector3(1./7., 1./11., 1./3.), /* test point */
        i, /* order of h1 field */
        i, /* order of test field */
        apf::VECTOR);
  }

  // test fields interpolating a user-defined matrix field
  for (int i = 1; i <= 6; i++) {
    if(0==PCUObj.Self())
      lion_oprint(1, "----TESTING MATRIX FIELD OF ORDER %d----\n", i);
    testH1(
        m, /* mesh */
        apf::Vector3(1./7., 1./11., 1./3.), /* test point */
        i, /* order of h1 field */
        i, /* order of test field */
        apf::MATRIX);
  }

  apf::destroyMesh(m);
  }
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

void M_exact(const apf::Vector3& x, apf::Matrix3x3& mat, int p)
{
  // Polynomial coefficients for each component of exact vector field
  double a[7] = { 1.0, -1.0,  2., -2., -1.0,  1.0, -1.0};
  double b[7] = {-2.0,  1.0, -2.,  2., -1.0, -1.0,  1.0};
  double c[7] = { 3.0,  0.0, -1.,  100., -1.0,  1.0,  0.0};

  apf::Vector3 value;
  value[0] = 0.0;
  value[1] = 0.0;
  value[2] = 0.0;
  for (int i = p;  i >= 0; i--) {
    value[0] += pow(x[0],p)*a[p];
    value[1] += pow(x[1],p)*b[p];
    value[2] += pow(x[2],p)*c[p];
  }

  mat[0] = value;
  mat[1] = apf::Vector3(value[2], 1., 0.);
  mat[2] = apf::Vector3(0., value[0], 1.);
}



void testH1(
    apf::Mesh2* m,
    const apf::Vector3& testXi,
    int ndOrder, int exactOrder,
    int dataType)
{
  PCU_ALWAYS_ASSERT(dataType == apf::VECTOR || dataType == apf::MATRIX);
  // Loop over all nodes and set scalar dofs.
  int dim = m->getDimension();
  apf::Field* h1Field;

  if (dataType == apf::VECTOR)
    h1Field = apf::createField(
      m, "h1_test", apf::VECTOR, apf::getH1Shape(ndOrder));
  else
    h1Field = apf::createField(
      m, "h1_test", apf::MATRIX, apf::getH1Shape(ndOrder));

  apf::MeshEntity* ent;
  apf::MeshIterator* it;

  for (int d = 0; d <= dim; d++) {
    if (!h1Field->getShape()->countNodesOn(apf::Mesh::simplexTypes[d])) {
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
      int non = h1Field->getShape()->countNodesOn(type);
      apf::MeshElement* me = apf::createMeshElement(m, ent);
      for (int i = 0; i < non; i++)
      {
        apf::Vector3 xi, p;
        h1Field->getShape()->getNodeXi(type, i, xi);
        apf::mapLocalToGlobal(me, xi, p);
        if (dataType == apf::VECTOR) {
          apf::Vector3 value;
          E_exact(p, value, exactOrder);
          apf::setVector(h1Field, ent, i, value);
        }
        else {
          apf::Matrix3x3 mat;
          M_exact(p, mat, exactOrder);
          apf::setMatrix(h1Field, ent, i, mat);
        }
      }
      apf::destroyMeshElement(me);
    }
    m->end(it);
  }


  // Verify that interpolated solution field agrees with exact field.
  for (int d = 0; d <= dim; d++) {

    double L2Error = 0.;
    it = m->begin(d);
    while( (ent = m->iterate(it)) ) {
      apf::MeshElement* me = apf::createMeshElement(m, ent);
      apf::Vector3 x;
      apf::mapLocalToGlobal(me, testXi, x);
      apf::Element* el = apf::createElement(h1Field, me);
      double err;
      if (dataType == apf::VECTOR) {
        apf::Vector3 eFieldExact;
        E_exact(x, eFieldExact, exactOrder);
        apf::Vector3 eFieldValue;
        apf::getVector(el, testXi, eFieldValue);
        err = ((eFieldValue - eFieldExact) * (eFieldValue - eFieldExact));
        err /= (eFieldExact * eFieldExact); // normalization factor
      }
      else {
        apf::Matrix3x3 mFieldExact;
        M_exact(x, mFieldExact, exactOrder);
        apf::Matrix3x3 mFieldValue;
        apf::getMatrix(el, testXi, mFieldValue);
        err = apf::getInnerProduct(
            mFieldValue - mFieldExact,
            mFieldValue - mFieldExact);
        err /= apf::getInnerProduct(mFieldExact, mFieldExact);
      }
      L2Error += err;
      apf::destroyElement(el);
      apf::destroyMeshElement(me);
    }
    m->end(it);

    // check for field interpolation
    if(0==m->getPCU()->Self()) {
      lion_oprint(1, "L2Error for entities of dimension %d is %e\n", d, L2Error);
    }
    PCU_ALWAYS_ASSERT_VERBOSE(L2Error < 1.e-12,
	"Fields were not interpolated correctly!");

  }
  m->removeField(h1Field);
  apf::destroyField(h1Field);
}
