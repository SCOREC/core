#include <dsp.h>
#include <apf.h>
#include <gmi_mesh.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <maSize.h>
#include <maMesh.h>
#include <PCU.h>
#include <sstream>
#include <vector>

using namespace std;

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
  fixed.insert(m->findModelEntity(2, 26));
  fixed.insert(m->findModelEntity(2, 6));
  fixed.insert(m->findModelEntity(2, 17));
  dsp::closeBoundary(m, fixed);
//  dsp::Smoother* smoother = dsp::Smoother::makeSemiSpring();
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
  vector < apf::MeshEntity* > V_total;
  int in_0; int fb_0;
  smoother->preprocess(m, fixed, moving, V_total, in_0, fb_0);
  /* number of displacement steps */
  for (int i = 0; i < 10; ++i) {
    apf::Field* dsp = dsp::applyRigidMotion(m, moving, r, t);
    smoother->smooth(dsp, V_total, in_0, fb_0);
  //dsp::tryToDisplace(m, dsp);
    apf::axpy(1, dsp, m->getCoordinateField());
    apf::destroyField(dsp);
    writeStep(m, i + 1);
    smoother->cleanup();
  }
  delete smoother;
  delete adapter;
  m->destroyNative();
  apf::destroyMesh(m);
  PCU_Comm_Free();
  MPI_Finalize();
}

