#include "ma.h"
#include "samSz.h"
#include "sam.h"
#include <chef.h>
#include <apf.h>
#include <gmi_mesh.h>
#include <apfMDS.h>
#include <apfShape.h>
#include <PCU.h>
#include <pcu_io.h>
#include <phRestart.h>
#include <phPartition.h>
#ifdef HAVE_SIMMETRIX
#include <gmi_sim.h>
#include <SimUtil.h>
#endif
#include <cassert>

static bool overwriteAPFCoord(apf::Mesh2* m) {
  apf::Field* f = m->findField("motion_coords");
  assert(f);
  double* vals = new double[apf::countComponents(f)];
  assert(apf::countComponents(f) == 3);
  apf::MeshEntity* vtx;
  apf::Vector3 points;
  apf::MeshIterator* itr = m->begin(0);
  while( (vtx = m->iterate(itr)) ) {
   apf::getComponents(f, vtx, 0, vals);
    for ( int i = 0; i < 3; i++ )  points[i] = vals[i];
    m->setPoint(vtx, 0, points);
  }
  m->end(itr);
  delete [] vals;
  return true;
}

static FILE* openfile_read(ph::Input&, const char* path) {
  return pcu_group_open(path, false);
}

int main(int argc, char** argv)
{
  assert(argc==3);
  const char* modelFile = argv[1];
  const char* meshFile = argv[2];
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
#ifdef HAVE_SIMMETRIX
  SimUtil_start();
  Sim_readLicenseFile(0);
  gmi_sim_start();
  gmi_register_sim();
#endif
  gmi_register_mesh();
  /* load model, mesh and configure input */
  apf::Mesh2* m = apf::loadMdsMesh(modelFile,meshFile);
  m->verify();
  ph::Input in;
  in.load("adapt.inp");
  in.openfile_read = openfile_read;
  /* attach solution and other fields */
  ph::readAndAttachFields(in,m);
  assert(overwriteAPFCoord(m));
  if (m->findField("material_type"))
    apf::destroyField(m->findField("material_type"));
  if (m->findField("meshQ"))
    apf::destroyField(m->findField("meshQ"));
  /* prepare size field */
  apf::Field* isoFld = samSz::isoSize(m);
  apf::Field* szFld = sam::multiplySF(m, isoFld, 1.0);
  apf::writeVtkFiles("before",m);
  /* mesh adaptation */
  ma::Input* ma_in = ma::configure(m, szFld);
  ma_in->shouldRunPreZoltan = true;
  ma_in->shouldRunMidParma = true;
  ma_in->shouldRunPostParma = true;
  ma_in->shouldSnap = in.snap;
  ma_in->shouldTransferParametric = in.transferParametric;
  if (m->hasMatching()) {
    ma_in->shouldSnap = false;
    ma_in->shouldFixShape = false;
  }
  ma::adapt(ma_in);
  m->verify();
  apf::writeVtkFiles("after",m);
  if (in.prePhastaBalanceMethod != "none" && PCU_Comm_Peers() > 1)
    ph::balance(in,m);
  /* output restart and geombc */
  chef::preprocess(m,in);
  m->destroyNative();
  apf::destroyMesh(m);
#ifdef HAVE_SIMMETRIX
  gmi_sim_stop();
  Sim_unregisterAllKeys();
  SimUtil_stop();
#endif
  PCU_Comm_Free();
  MPI_Finalize();
}

