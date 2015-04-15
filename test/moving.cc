#include <dsp.h>
#include <apf.h>
#include <gmi_mesh.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <maSize.h>
#include <maMesh.h>
#include <PCU.h>
#include <sstream>

static void writeStep(apf::Mesh* m, int i)
{
  std::stringstream ss;
  ss << "step_" << i << "_";
  std::string s = ss.str();
  apf::writeVtkFiles(s.c_str(), m);
}

int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  if ( argc != 3 ) {
    fprintf(stderr, "Usage: %s <model> <mesh>\n", argv[0]);
    return 0;
  }
  gmi_register_mesh();
  apf::Mesh2* m = apf::loadMdsMesh(argv[1],argv[2]);
  dsp::Boundary moving;
  moving.insert(m->findModelEntity(2, 57));
  moving.insert(m->findModelEntity(2, 62));
  moving.insert(m->findModelEntity(2, 66));
  dsp::closeBoundary(m, moving);
  dsp::Boundary fixed;
  dsp::Smoother* smoother = dsp::Smoother::makeLaplacian();
//double avgEdgeLen = ma::getAverageEdgeLength(m);
//dsp::Adapter* adapter = dsp::Adapter::makeUniform(avgEdgeLen);
  dsp::Adapter* adapter = dsp::Adapter::makeEmpty();
  apf::Vector3 a;
  apf::Vector3 b;
  ma::getBoundingBox(m, a, b);
  double zDist = b.z() - a.z();
  apf::Vector3 t(0,0, zDist / 30);
  apf::Matrix3x3 r(1,0,0,
                   0,1,0,
                   0,0,1);
  writeStep(m, 0);
  /* number of displacement steps */
  for (int i = 0; i < 1; ++i) {
    apf::Field* dsp = dsp::applyRigidMotion(m, moving, r, t);
    smoother->smooth(dsp, fixed, moving);
  //dsp::tryToDisplace(m, dsp);
    apf::axpy(1, dsp, m->getCoordinateField());
    apf::destroyField(dsp);
    writeStep(m, i + 1);
  }
  delete smoother;
  delete adapter;
  m->destroyNative();
  apf::destroyMesh(m);
  PCU_Comm_Free();
  MPI_Finalize();
}

