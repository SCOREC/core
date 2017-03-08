#include <PCU.h>
#include <apf.h>
#include <apfMDS.h>
#include <apfBox.h>
#include <apfShape.h>
#include <apfNumbering.h>
#include <mthQR.h>
#include <pcu_util.h>
#include <cstdlib>
#include <iomanip>
#include <fstream>

namespace {

static double const pi = 3.14159265358979323846;

static double forcing(apf::Vector3 const& p, int d) {
  double x = p[0];
  double y = p[1];
  double z = p[2];
  double v = 4.0*d*pi*pi*std::sin(2.0*pi*x);
  if (d > 1) v *= std::sin(2.0*pi*y);
  if (d > 2) v *= std::sin(2.0*pi*z);
  return v;
}

static double solution(apf::Vector3 const& p, int d) {
  double x = p[0];
  double y = p[1];
  double z = p[2];
  double v = std::sin(2.0*pi*x);
  if (d > 1) v *= std::sin(2.0*pi*y);
  if (d > 2) v *= std::sin(2.0*pi*z);
  return v;
}

class Poisson {

  public:

    Poisson(int d, int p, int n) {
      num_dims = d;
      p_order = p;
      q_degree = 6;
      num_grid = n;
      x = apf::Vector3(0,0,0);
      xi = apf::Vector3(0,0,0);
    }

    ~Poisson() {
      apf::destroyField(sol);
      apf::destroyField(sol_mid);
      apf::destroyField(sol_vtx);
      apf::destroyGlobalNumbering(numbering);
      mesh->destroyNative();
      apf::destroyMesh(mesh);
    }

    void run() {
      setup_grid();
      setup_lin_alg();
      fill_volumetric();
      fill_boundary();
      solve_lin_sys();
      attach_nodal_solution();
      attach_mid_solution();
      attach_vtx_solution();
      compute_l2_error();
    }

  private:

    int num_dims;
    int p_order;
    int q_degree;
    int num_grid;
    int num_dofs;

    apf::Mesh2* mesh;
    apf::Field* sol;
    apf::Field* sol_mid;
    apf::Field* sol_vtx;
    apf::FieldShape* shape;
    apf::GlobalNumbering* numbering;

    mth::Matrix<double> K;
    mth::Vector<double> u;
    mth::Vector<double> F;

    apf::Vector3 x;
    apf::Vector3 xi;

    apf::NewArray<double> BF;
    apf::NewArray<apf::Vector3> GBF;
    apf::NewArray<long> numbers;

    void setup_grid() {
      int nx = num_grid;
      int ny = (num_dims > 1) ? num_grid : 0;
      int nz = (num_dims > 2) ? num_grid : 0;
      double wx = 1.0;
      double wy = (num_dims > 1) ? 1.0 : 0.0;
      double wz = (num_dims > 2) ? 1.0 : 0.0;
      mesh = apf::makeMdsBox(nx, ny, nz, wx, wy, wz, true);
      apf::reorderMdsMesh(mesh);
      shape = apf::getLagrange(p_order);
      sol = apf::createField(mesh, "u", apf::SCALAR, shape);
      sol_mid = apf::createStepField(mesh, "u_mid", apf::SCALAR);
      sol_vtx = apf::createFieldOn(mesh, "u_vtx", apf::SCALAR);
      numbering = apf::makeGlobal(apf::numberOwnedNodes(mesh, "n", shape));
    }

    void setup_lin_alg() {
      num_dofs = apf::countNodes(numbering);
      K.resize(num_dofs, num_dofs);
      x.resize(num_dofs);
      F.resize(num_dofs);
      K.zero();
      K.zero();
      F.zero();
    }

    void fill_volumetric() {
      apf::MeshEntity* elem;
      apf::MeshIterator* elems = mesh->begin(num_dims);
      while ((elem = mesh->iterate(elems))) {
        apf::MeshElement* me = apf::createMeshElement(mesh, elem);
        int num_ips = apf::countIntPoints(me, q_degree);
        for (int ip=0; ip < num_ips; ++ip) {
          apf::getIntPoint(me, q_degree, ip, xi);
          apf::mapLocalToGlobal(me, xi, x);
          apf::getBF(shape, me, xi, BF);
          apf::getGradBF(shape, me, xi, GBF);
          int nn = apf::getElementNumbers(numbering, elem, numbers);
          double w = apf::getIntWeight(me, q_degree, ip);
          double dv = apf::getDV(me, xi);
          double f = forcing(x, num_dims);
          for (int i=0; i < nn; ++i)
            F( numbers[i] ) += f * BF[i] * w * dv;
          for (int i=0; i < nn; ++i)
          for (int j=0; j < nn; ++j)
          for (int d=0; d < num_dims; ++d)
            K( numbers[i], numbers[j] ) += GBF[j][d] * GBF[i][d] * w * dv;
        }
        apf::destroyMeshElement(me);
      }
      mesh->end(elems);
    }

    void fill_boundary() {
      gmi_model* model = mesh->getModel();
      gmi_ent* boundary;
      gmi_iter* boundaries = gmi_begin(model, num_dims-1);
      while ((boundary = gmi_next(model, boundaries))) {
        apf::DynamicArray<apf::Node> nodes;
        apf::ModelEntity* bdry =
          reinterpret_cast<apf::ModelEntity*>(boundary);
        apf::getNodesOnClosure(mesh, bdry, nodes, shape);
        for (size_t n=0; n < nodes.getSize(); ++n) {
          int row = apf::getNumber(numbering, nodes[n]);
          F(row) = 0.0;
          for (int col=0; col < num_dofs; ++col)
            K(row, col) = 0.0;
          K(row, row) = 1.0;
        }
      }
      gmi_end(model, boundaries);
    }

    void solve_lin_sys() {
      bool solved = mth::solveQR(K, F, u);
      if (! solved) apf::fail("QR failed!");
    }

    void attach_nodal_solution() {
      apf::DynamicArray<apf::Node> nodes;
      apf::getNodes(numbering, nodes);
      for (size_t i=0; i < nodes.size(); ++i) {
        int n = apf::getNumber(numbering, nodes[i]);
        apf::setScalar(sol, nodes[i].entity, nodes[i].node, u(n));
      }
    }

    void attach_mid_solution() {
      apf::MeshEntity* elem;
      apf::MeshIterator* elems = mesh->begin(num_dims);
      if (num_dims == 1) xi = apf::Vector3(0, 0, 0);
      if (num_dims == 2) xi = apf::Vector3(1.0/3.0, 1.0/3.0, 0);
      if (num_dims == 3) xi = apf::Vector3(0.25, 0.25, 0);
      while ((elem = mesh->iterate(elems))) {
        apf::MeshElement* me = apf::createMeshElement(mesh, elem);
        apf::Element* e = apf::createElement(sol, me);
        double u = apf::getScalar(e, xi);
        apf::setScalar(sol_mid, elem, 0, u);
        apf::destroyElement(e);
        apf::destroyMeshElement(me);
      }
      mesh->end(elems);
    }

    void attach_vtx_solution() {
      apf::MeshEntity* vtx;
      apf::MeshIterator* vertices = mesh->begin(0);
      while ((vtx = mesh->iterate(vertices))) {
        double uh = apf::getScalar(sol, vtx, 0);
        apf::setScalar(sol_vtx, vtx, 0, uh);
      }
      mesh->end(vertices);
    }

    void compute_l2_error() {
      int high_q = 6;
      double L2_error = 0.0;
      apf::MeshEntity* elem;
      apf::MeshIterator* elems = mesh->begin(num_dims);
      while ((elem = mesh->iterate(elems))) {
        apf::MeshElement* me = apf::createMeshElement(mesh, elem);
        apf::Element* e = apf::createElement(sol, me);
        int num_ips = apf::countIntPoints(me, high_q);
        for (int ip=0; ip < num_ips; ++ip) {
          apf::getIntPoint(me, 6, ip, xi);
          apf::mapLocalToGlobal(me, xi, x);
          double w = apf::getIntWeight(me, high_q, ip);
          double dv = apf::getDV(me, xi);
          double u = solution(x, num_dims);
          double uh = apf::getScalar(e, xi);
          L2_error += std::abs(u-uh) * std::abs(u-uh) * w * dv;
        }
        apf::destroyElement(e);
        apf::destroyMeshElement(me);
      }
      L2_error = std::sqrt(L2_error);
      printf("n grid: %d | n dofs: %d | L2 error: %.15e\n",
          num_grid, num_dofs, L2_error);
    }

};

void test(int dim, int p) {
  int steps = 6-dim;
  int n_grid = 6-dim;
  for (int i=0; i < steps; ++i) {
    Poisson poisson(dim, p, n_grid);
    poisson.run();
    n_grid *= 2;
  }
}

}

int main(int argc, char** argv) {
  PCU_ALWAYS_ASSERT(argc == 3);
  MPI_Init(&argc, &argv);
  PCU_Comm_Init();
  PCU_ALWAYS_ASSERT(! PCU_Comm_Self());
  int dim = atoi(argv[1]);
  int p = atoi(argv[2]);
  test(dim, p);
  PCU_Comm_Free();
  MPI_Finalize();
}
