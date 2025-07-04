#include <iostream>
#include <iomanip>

#include <PCU.h>
#include <apf.h>
#include <apfCAP.h>
#include <apfMesh2.h>
#include <gmi_cap.h>
#include <lionPrint.h>
#include <pcu_util.h>

void checkParametrization(apf::Mesh* mesh);

int main(int argc, char** argv)
{
  pcu::Init(&argc, &argv);
  { // pcu object scope
  pcu::PCU PCUobj;

  if (argc != 2) {
    if (PCUobj.Self() == 0)
      std::cerr << "usage: " << argv[0] << " <cre file .cre>\n";
    return EXIT_FAILURE;
  }

  const char* creFileName = argv[1];

  lion_set_verbosity(1);
  gmi_cap_start();
  gmi_register_cap();

  gmi_model* model = gmi_cap_load(creFileName);
  apf::Mesh2* mesh = apf::createCapMesh(model, &PCUobj);

  // check parametrization using capstone apis
  checkParametrization(mesh);
  apf::destroyMesh(mesh);
  gmi_cap_stop();
  } // pcu object scope
  pcu::Finalize();
}

void checkParametrization(apf::Mesh* mesh) {
  int count[2] = {0, 0};
  double sum = 0.0;
  std::cout <<
    "ct, (        u,         v), (     umin,      umax), "
    "(     vmin,      vmax), diff\n";
  apf::MeshIterator* it = mesh->begin(0);
  for (apf::MeshEntity* vtx; (vtx = mesh->iterate(it));) {
    apf::ModelEntity* me = mesh->toModel(vtx);
    int dim = mesh->getModelType(me);
    if (dim != 1 && dim != 2) continue;
    double range_u[2], range_v[2];
    mesh->getPeriodicRange(me, 0, range_u);
    mesh->getPeriodicRange(me, 1, range_v);
    apf::Vector3 xi;
    mesh->getParam(vtx, xi);
    // coordinate from mesh
    apf::Vector3 coord;
    mesh->getPoint(vtx, 0, coord);
    // coordinate from surface
    apf::Vector3 pcoord;
    mesh->snapToModel(me, xi, pcoord);
    apf::Vector3 diff = coord - pcoord;

    constexpr int dim_printmax = 25;
    if (count[dim - 1] < dim_printmax) {
      std::cout << std::setw(2) << count[dim - 1] << ", "
        << std::setprecision(3) << std::scientific
        << '(' << std::setw(8) << xi.x() << ", "
        << std::setw(8) << xi.y() << "), "
        << '(' << std::setw(4) << range_u[0] << ", "
        << std::setw(4) << range_u[1] << "), ";
      if (dim == 2)
        std::cout << '(' << std::setw(4) << range_v[0] << ", "
          << std::setw(4) << range_v[1] << "), ";
      else
        std::cout << "----------------------, ";
      std::cout << std::setw(4) << std::defaultfloat << diff.getLength()
        << std::endl;
    } else if (count[dim - 1] == dim_printmax) {
      std::cout << "Skipping printing remaining "
        << (dim == 1 ? "edges" : "faces") << std::endl;
    }
    sum += diff * diff;
    count[dim - 1]++;
  }
  std::cout << "norm of the difference vector is " << std::sqrt(sum)
    << std::endl;
}
