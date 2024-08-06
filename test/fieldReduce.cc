#include <gmi_mesh.h>
#include <apf.h>
#include <apfMesh2.h>
#include <apfMDS.h>
#include <apfShape.h>
#include <lionPrint.h>
#include <parma.h>
#ifdef HAVE_SIMMETRIX
#include <gmi_sim.h>
#include <SimUtil.h>
#include <MeshSim.h>
#include <SimModel.h>
#endif
#include <pcu_util.h>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <memory>

namespace {

const char* modelFile = 0;
const char* meshFile = 0;
const char* outFile = 0;
int myrank = -1;
int commsize = -1;

double getValue(apf::Vector3& coords, double addval)
{
  return coords.x() + coords.y() + 1 + addval;
}

apf::Field* getTestField(apf::Mesh* m, const char* fname, double addval)
{
  apf::FieldShape* fshape = apf::getLagrange(1);
  apf::Field* f = apf::createLagrangeField(m, fname, apf::SCALAR, 1);

  // write known data to field, check that the sum is n times the number of
  // copies
  apf::MeshEntity* e;
  apf::MeshIterator* it;
  apf::Vector3 coords;

  for (int d=0; d <= m->getDimension(); ++d)
  {
    if (!fshape->hasNodesIn(d))
      continue;

    it = m->begin(d);
    while ( (e = m->iterate(it)) )
    {
      m->getPoint(e, 0, coords);
      double val = getValue(coords, addval);
      apf::setScalar(f, e, 0, val);
    }  // end while
    m->end(it);
  }  // end for

  return f;
}

bool testReduce(apf::Mesh* m, int casenum)
{

  double addval = 0;
  if (casenum == 0)  // sum
    addval = 0;
  else
    addval = myrank;

  char fname[256];
  snprintf(fname, 256, "test%d", casenum);
  apf::Field* f = getTestField(m, fname, addval);
  apf::FieldShape* fshape = apf::getShape(f);
  apf::Sharing* shr = apf::getSharing(m);

  if (casenum == 0)
    apf::accumulate(f, shr);
  else if (casenum == 1)
    apf::sharedReduction(f, shr, false, apf::ReductionMax<double>());
  else if (casenum == 2)
    apf::sharedReduction(f, shr, false, apf::ReductionMin<double>());

  // verify the result is n times the number of copies
  apf::MeshEntity* e;
  apf::MeshIterator* it;
  apf::Vector3 coords;
  bool failflag = false;

  for (int d=0; d <= m->getDimension(); ++d)
  {
    if (!fshape->hasNodesIn(d))
      continue;

    it = m->begin(d);
    while ( (e = m->iterate(it)) )
    {
      apf::CopyArray copies;  // need to redeclare this every iteration so
                              // the size gets reset
      shr->getCopies(e, copies);
      int ntimes = copies.getSize() + 1;

      if (ntimes > 1)
      {
        // this entity has remotes, get some info
        m->getPoint(e, 0, coords);
        double val_f = apf::getScalar(f, e, 0);

        // get max and min peer
        int maxpeer = myrank;
        int minpeer = myrank;
        for (std::size_t j = 0; j < copies.getSize(); ++j)
        {
          if (copies[j].peer > maxpeer)
            maxpeer = copies[j].peer;

          if (copies[j].peer < minpeer)
            minpeer = copies[j].peer;
        }

        // do the test
        if (casenum == 0)
        {
          double val = getValue(coords, addval);
          failflag = ( failflag || ( std::fabs(ntimes*val - val_f) > 1e-13 ) );
        } else if (casenum == 1) // max
        {
          double val = getValue(coords, maxpeer);
          failflag = ( failflag || ( std::fabs(val - val_f) > 1e-13 ) );
        } else if (casenum == 2) { // min
          double val = getValue(coords, minpeer);
          failflag = ( failflag || ( std::fabs(val - val_f) > 1e-13 ) );
        }  // end if casenum

      }  // end if ntimes
    } // end while
    m->end(it);
  }  // end for

  delete shr;
  return failflag;
}

void freeMesh(apf::Mesh* m)
{
  m->destroyNative();
  apf::destroyMesh(m);
}

void getConfig(int argc, char** argv, pcu::PCU *PCUObj)
{
  if ( argc != 3 ) {
    if ( !PCUObj->Self() )
      printf("Usage: %s <model> <mesh>\n", argv[0]);
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }
  modelFile = argv[1];
  meshFile = argv[2];
  outFile = argv[3];
//  PCU_ALWAYS_ASSERT(commsize == 4);  //DEBUGGING
}

}

int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  bool failflag = false;
  {
  auto PCUObj = std::unique_ptr<pcu::PCU>(new pcu::PCU(MPI_COMM_WORLD));
  lion_set_verbosity(1);
  MPI_Comm_rank(PCUObj.get()->GetMPIComm(), &myrank);
  MPI_Comm_size(PCUObj.get()->GetMPIComm(), &commsize);
#ifdef HAVE_SIMMETRIX
  MS_init();
  SimModel_start();
  Sim_readLicenseFile(0);
  gmi_sim_start();
  gmi_register_sim();
#endif
  gmi_register_mesh();
  getConfig(argc,argv,PCUObj.get());
  gmi_model* g = 0;
  g = gmi_load(modelFile);
  apf::Mesh2* m = 0;
  m = apf::loadMdsMesh(g, meshFile, PCUObj.get());

  
  for (int i=0; i < 3; ++i)
    failflag = failflag || testReduce(m, i);

  freeMesh(m);
#ifdef HAVE_SIMMETRIX
  gmi_sim_stop();
  Sim_unregisterAllKeys();
  SimModel_stop();
  MS_exit();
#endif
  }
  MPI_Finalize();

  return failflag;
  
}

