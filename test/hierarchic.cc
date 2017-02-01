#include <PCU.h>
#include <apf.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <apfShape.h>
#include <apfNumbering.h>
#include <gmi_mesh.h>
#include <mthQR.h>
#include <cassert>
#include <cstdlib>

#include <iostream>
#include <iomanip>
#include <fstream>

namespace {

static double linear(int d, apf::Vector3 const& p) {
  double x = p[0];
  double y = p[1];
  double z = p[2];
  double v = 1.0 + x;
  if (d > 1) v += y;
  if (d > 2) v += z;
  return v;
}

static double quadratic(int d, apf::Vector3 const& p) {
  double x = p[0];
  double y = p[1];
  double z = p[2];
  double v = 1.0 + x + x*x;
  if (d > 1) v += y + x*y + y*y;
  if (d > 2) v += z + x*z + y*z + z*z;
  return v;
}

static double cubic(int d, apf::Vector3 const& p) {
  double x = p[0];
  double y = p[1];
  double z = p[2];
  double v = 1.0 + x*x*x;
  if (d > 1) v += x*y*y + y*y*y;
  if (d > 2) v += z*z*z;
  return v;
}

static double function(int p_order, int d, apf::Vector3 const& p) {
  if (p_order == 1) return linear(d, p);
  else if (p_order == 2) return quadratic(d, p);
  else if (p_order == 3) return cubic(d, p);
  else apf::fail("invalid polynomial order\n");
}

class L2Projector {
  public:
    L2Projector(apf::Mesh* m, int p);
    void run();
  private:
    void initialize();
    void fill();
    void write_matrix();
    void write_vector();
    void solve();
    void attach();
    void compare();
    void finalize();
    apf::Mesh* mesh;
    apf::Field* field;
    apf::FieldShape* shape;
    apf::GlobalNumbering* numbering;
    mth::Matrix<double> M;
    mth::Vector<double> x;
    mth::Vector<double> b;
    int p_order;
    int q_degree;
};

L2Projector::L2Projector(apf::Mesh* m, int p) {
  mesh = m;
  field = 0;
  shape = 0;
  numbering = 0;
  p_order = p;
  q_degree = 2 * p_order;
  printf("L2 projector: polynomial order: %d\n", p_order);
  printf("L2 projector: quadrature degree: %d\n", q_degree);
}

void L2Projector::initialize() {
  printf("L2 projector: initializing\n");
  shape = apf::getHierarchic(p_order);
  field = apf::createField(mesh, "u", apf::SCALAR, shape);
  numbering = apf::makeGlobal(apf::numberOwnedNodes(mesh, "n", shape));
  int N = countNodes(numbering);
  M.resize(N,N);
  x.resize(N);
  b.resize(N);
  M.zero();
  x.zero();
  b.zero();
  printf("L2 projector: number dofs: %d\n", N);
}

void L2Projector::fill() {
  printf("L2 projector: fill\n");
  apf::Vector3 x(0,0,0);
  apf::Vector3 xi(0,0,0);
  apf::NewArray<long> numbers;
  apf::NewArray<double> BF;
  apf::MeshEntity* elem;
  int d = mesh->getDimension();
  apf::MeshIterator* elems = mesh->begin(d);
  while ((elem = mesh->iterate(elems))) {
    apf::MeshElement* me = apf::createMeshElement(mesh, elem);
    int num_ips = apf::countIntPoints(me, q_degree);
    for (int ip=0; ip < num_ips; ++ip) {
      apf::getIntPoint(me, q_degree, ip, xi);
      apf::mapLocalToGlobal(me, xi, x);
      apf::getBF(shape, me, xi, BF);
      int num_nodes = apf::getElementNumbers(numbering, elem, numbers);
      double w = apf::getIntWeight(me, q_degree, ip);
      double dv = apf::getDV(me, xi);
      double u = function(p_order, d, x);
      for (int i=0; i < num_nodes; ++i)
        b( numbers[i] ) += u * BF[i] * w * dv;
      for (int i=0; i < num_nodes; ++i)
      for (int j=0; j < num_nodes; ++j)
        M( numbers[i], numbers[j] ) += BF[i] * BF[j] * w * dv;
    }
    apf::destroyMeshElement(me);
  }
  mesh->end(elems);
}

void L2Projector::solve() {
  printf("L2 projector: solving\n");
  bool solved = mth::solveQR(M, b, x);
  if (! solved) apf::fail("QR failed!");
}

void L2Projector::attach() {
  printf("L2 projector: attaching to field\n");
  apf::DynamicArray<apf::Node> nodes;
  apf::getNodes(numbering, nodes);
  for (size_t i=0; i < nodes.size(); ++i) {
    apf::MeshEntity* ent = nodes[i].entity;
    int n = nodes[i].node;
    long N = apf::getNumber(numbering, nodes[i]);
    apf::setScalar(field, ent, n, x(N));
  }
}

void L2Projector::compare() {
  printf("L2 projector: compare\n");
  apf::MeshEntity* elem;
  apf::Vector3 x(0,0,0);
  apf::Vector3 p(0.25,0,0);
  int d = mesh->getDimension();
  if (d > 1) p[1] = 0.25;
  if (d > 2) p[2] = 0.25;
  apf::MeshIterator* elems = mesh->begin(d);
  while ((elem = mesh->iterate(elems))) {
    apf::MeshElement* me = apf::createMeshElement(mesh, elem);
    apf::Element* e = apf::createElement(field, me);
    apf::mapLocalToGlobal(me, p, x);
    double v = function(p_order, d, x);
    double vh = apf::getScalar(e, p);
    assert(std::abs(v - vh) < 1.0e-13);
    apf::destroyElement(e);
    apf::destroyMeshElement(me);
  }
  mesh->end(elems);
}

void L2Projector::finalize() {
  apf::destroyField(field);
  apf::destroyGlobalNumbering(numbering);
}

void L2Projector::run() {
  initialize();
  fill();
  solve();
  attach();
  compare();
  finalize();
}

void test(apf::Mesh* m, int p_order) {
  L2Projector projector(m, p_order);
  projector.run();
  printf("test passed: that's real neat-o!\n");
}

}

int main(int argc, char** argv)
{
  assert(argc==4);
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  assert(! PCU_Comm_Self());
  gmi_register_mesh();
  apf::Mesh2* m = apf::loadMdsMesh(argv[1], argv[2]);
  apf::reorderMdsMesh(m);
  m->verify();
  int p_order = atoi(argv[3]);
  test(m, p_order);
  m->destroyNative();
  apf::destroyMesh(m);
  PCU_Comm_Free();
  MPI_Finalize();
}
