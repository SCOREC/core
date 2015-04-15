#include <dsp.h>
#include <apf.h>
#include <gmi_mesh.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <maSize.h>
#include <maMesh.h>
#include <PCU.h>

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
  dsp::Boundary fixed;
//dsp::Smoother* smoother = dsp::Smoother::makeLagrangian();
  dsp::Smoother* smoother = dsp::Smoother::makeEmpty();
//double avgEdgeLen = ma::getAverageEdgeLength(m);
//dsp::Adapter* adapter = dsp::Adapter::makeUniform(avgEdgeLen);
  dsp::Adapter* adapter = dsp::Adapter::makeEmpty();
  apf::Vector3 a;
  apf::Vector3 b;
  ma::getBoundingBox(m, a, b);
  double zDist = b.z() - a.z();
  apf::Vector3 t(0,0, zDist / 10);
  apf::Matrix3x3 r(1,0,0,
                   0,1,0,
                   0,0,1);
  /* number of displacement steps */
  for (int i = 0; i < 1; ++i) {
    apf::Field* dsp = dsp::applyRigidMotion(m, fixed, r, t);
    dsp::displace(m, dsp, smoother, adapter, fixed, moving);
    apf::destroyField(dsp);
  }
  delete smoother;
  delete adapter;
  apf::writeVtkFiles("after", m);
  m->destroyNative();
  apf::destroyMesh(m);
  PCU_Comm_Free();
  MPI_Finalize();
}

