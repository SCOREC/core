#include "ma.h"
#include <apf.h>
#include <gmi_mesh.h>
#include <apfMDS.h>
#include <apfShape.h>
#include <PCU.h>

int main(int argc, char** argv)
{
  assert(argc==3);
  const char* modelFile = argv[1];
  const char* meshFile = argv[2];
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  gmi_register_mesh();
  ma::Mesh* m = apf::loadMdsMesh(modelFile,meshFile);
  m->verify();
  apf::Field* scales;
  scales = apf::createLagrangeField(m, "proteus_size_scale", apf::VECTOR, 1);
  apf::Field* frames;
  frames = apf::createLagrangeField(m, "proteus_size_frame", apf::MATRIX, 1);
  ma::Input* in = ma::configure(m, scales, frames);
  in->shouldRunPreParma = true;
  in->shouldRunMidParma = true;
  in->shouldRunPostParma = true;
  in->maximumIterations = 1;
  in->shouldFixShape = true;
  ma::adapt(in);
  m->verify();
  apf::writeVtkFiles("after",m);
  m->destroyNative();
  apf::destroyMesh(m);
  PCU_Comm_Free();
  MPI_Finalize();
}

