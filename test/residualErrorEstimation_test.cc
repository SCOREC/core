#include "ma.h"
#include <ree.h>
#include <apf.h>
#include <gmi_mesh.h>
#include <apfMDS.h>
#include <apfShape.h>
#include <lionPrint.h>
#ifdef HAVE_SIMMETRIX
#include <gmi_sim.h>
#include <SimUtil.h>
#include <MeshSim.h>
#include <SimModel.h>
#endif
#include <pcu_util.h>

void E_exact(const apf::Vector3 &x, apf::Vector3& E);
double computeElementExactError(apf::Mesh* mesh, apf::MeshEntity* e,
  apf::Field* f);

double freq = 1.0, kappa;
int dim = 3;

int main(int argc, char** argv)
{
  PCU_ALWAYS_ASSERT(argc==3);
  const char* modelFile = argv[1];
  const char* meshFile = argv[2];
  MPI_Init(&argc,&argv);
  {
  pcu::PCU PCUObj = pcu::PCU(MPI_COMM_WORLD);
  lion_set_verbosity(1);
#ifdef HAVE_SIMMETRIX
  MS_init();
  SimModel_start();
  Sim_readLicenseFile(0);
  gmi_sim_start();
  gmi_register_sim();
#endif
  gmi_register_mesh();
  ma::Mesh* m = apf::loadMdsMesh(modelFile,meshFile,&PCUObj);
  m->verify();

  kappa = freq * M_PI;

  // get electric field
  apf::Field* electric_field = m->getField(0);
  PCU_ALWAYS_ASSERT(electric_field);

  // compute residual error
  apf::Field* residual_error_field = ree::estimateError(electric_field);

  // compute exact error
  apf::Field* exact_error_field = apf::createIPField(
      m, "exact_error_field", apf::SCALAR, 1);
  apf::MeshEntity* ent;
  apf::MeshIterator* itr = m->begin(3);
  while ((ent = m->iterate(itr)))
  {
    double exact_element_error = computeElementExactError(
        m, ent, electric_field);
    apf::setScalar(exact_error_field, ent, 0, exact_element_error);
  }
  m->end(itr);

  // get max and avg computed and exact errors
  double max_exact_error = 0., max_computed_error = 0.;
  double avg_exact_error = 0., avg_computed_error = 0.;
  itr = m->begin(3);
  while ((ent = m->iterate(itr)))
  {
    double exact_error = apf::getScalar(exact_error_field, ent, 0);
    if (exact_error > max_exact_error) max_exact_error = exact_error;
    avg_exact_error += exact_error;

    double computed_error = apf::getScalar(residual_error_field, ent, 0);
    if (computed_error > max_computed_error) max_computed_error = computed_error;
    avg_computed_error += computed_error;
  }
  m->end(itr);
  avg_exact_error /= m->count(3);
  avg_computed_error /= m->count(3);

  lion_oprint(1, "Max     Exact Error: %e\n", max_exact_error);
  lion_oprint(1, "Average Exact Error: %e\n", avg_exact_error);
  lion_oprint(1, "Max     Computed Error: %e\n", max_computed_error);
  lion_oprint(1, "Average Computed Error: %e\n", avg_computed_error);

  apf::destroyField(residual_error_field);
  apf::destroyField(exact_error_field);

  m->destroyNative();
  apf::destroyMesh(m);
#ifdef HAVE_SIMMETRIX
  gmi_sim_stop();
  Sim_unregisterAllKeys();
  SimModel_stop();
  MS_exit();
#endif
  }
  MPI_Finalize();
}

void E_exact(const apf::Vector3 &x, apf::Vector3& E)
{
   if (dim == 3)
   {
      E[0] = sin(kappa * x[1]);
      E[1] = sin(kappa * x[2]);
      E[2] = sin(kappa * x[0]);
   }
   else
   {
      E[0] = sin(kappa * x[1]);
      E[1] = sin(kappa * x[0]);
      E[2] = 0.0;
   }
}

double computeElementExactError(apf::Mesh* mesh, apf::MeshEntity* e,
  apf::Field* f)
{
  double error = 0.0;

  apf::FieldShape* fs = f->getShape();
  int type = mesh->getType(e);
  PCU_ALWAYS_ASSERT(type == apf::Mesh::TET);
  int dim = apf::getDimension(mesh, e);
  double w;

  apf::MeshElement* me = apf::createMeshElement(mesh, e);
  apf::Element* el = apf::createElement(f, me);
  int order = 2*fs->getOrder() + 1;
  int np = apf::countIntPoints(me, order);

  apf::Vector3 femsol, exsol;

  apf::Vector3 p;
  for (int i = 0; i < np; i++) {
    apf::getIntPoint(me, order, i, p);
    double weight = apf::getIntWeight(me, order, i);
    apf::Matrix3x3 J;
    apf::getJacobian(me, p, J);
    double jdet = apf::getJacobianDeterminant(J, dim);
    w = weight * jdet;

    apf::getVector(el, p, femsol);

    apf::Vector3 global;
    apf::mapLocalToGlobal(me, p, global);
    E_exact(global, exsol);
    apf::Vector3 diff = exsol - femsol;

    error += w * (diff * diff);
  }
  if (error < 0.0)
    error = -error;

  apf::destroyElement(el);
  apf::destroyMeshElement(me);

  return sqrt(error);
}

