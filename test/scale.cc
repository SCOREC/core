#include <apf.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <apfShape.h>
#include <gmi_mesh.h>
#include <lionPrint.h>
#include <cstdlib>

namespace {

struct Scale {
  double x;
  double y;
  double z;
  double s;
};

static void print_usage(char** argv, pcu::PCU *PCUObj) {
  if (! PCUObj->Self())
    printf("Usage: %s <model> <mesh> <out> <x> <y> <z> <scale>\n", argv[0]);
  MPI_Finalize();
  exit(EXIT_FAILURE);
}

static void scale_dim(
    int d, apf::Mesh2* m, apf::FieldShape* shape, Scale const& scale) {
  apf::Vector3 X(0,0,0);
  apf::MeshEntity* ent;
  apf::MeshIterator* it = m->begin(d);
  apf::Field* coords = m->getCoordinateField();
  while ((ent = m->iterate(it))) {
    int type = m->getType(ent);
    int nnodes = shape->countNodesOn(type);
    for (int n = 0; n < nnodes; ++n) {
      m->getPoint(ent, n, X);
      X[0] = (X[0] + scale.x) * scale.s;
      X[1] = (X[1] + scale.y) * scale.s;
      X[2] = (X[2] + scale.z) * scale.s;
      apf::setVector(coords, ent, n, X);
    }
  }
  m->end(it);
}

static void scale_mesh(apf::Mesh2* m, Scale const& s) {
  int dim = m->getDimension();
  apf::FieldShape* shape = m->getShape();
  for (int d = 0; d < dim; ++d)
    if (shape->hasNodesIn(d))
      scale_dim(d, m, shape, s);
}

}

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  {
  pcu::PCU PCUObj = pcu::PCU(MPI_COMM_WORLD);
  lion_set_verbosity(1);
  gmi_register_mesh();
  if (argc != 8) print_usage(argv, &PCUObj);
  const char* gfile = argv[1];
  const char* mfile = argv[2];
  const char* ofile = argv[3];
  Scale scale;
  scale.x = atof(argv[4]);
  scale.y = atof(argv[5]);
  scale.z = atof(argv[6]);
  scale.s = atof(argv[7]);
  apf::Mesh2* m = apf::loadMdsMesh(gfile, mfile, &PCUObj);
  m->verify();
  scale_mesh(m, scale);
  m->verify();
  m->writeNative(ofile);
  m->destroyNative();
  apf::destroyMesh(m);
  }
  MPI_Finalize();
}
