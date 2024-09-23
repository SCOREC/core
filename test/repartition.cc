#include <gmi_mesh.h>
#include <apf.h>
#include <apfMesh2.h>
#include <apfMDS.h>
#include <PCU.h>
#include <lionPrint.h>
#include <parma.h>
#include <apfZoltan.h>
#include <apfPartition.h>
#include <pcu_util.h>
#include <cstdlib>

namespace {

const char* modelFile = 0;
const char* meshFile = 0;
const char* outFile = 0;
int inputPartCount = 1;

struct CreateGroupCommResult{
  bool isOriginal;
  pcu::PCU *group_pcu_obj;
};

void freeMesh(apf::Mesh* m)
{
  m->destroyNative();
  apf::destroyMesh(m);
}

CreateGroupCommResult createGroupComm(pcu::PCU *PCUObj)
{
  apf::Contract contract(inputPartCount, PCUObj->Peers());
  int self = PCUObj->Self();
  PCU_ALWAYS_ASSERT(self == PCUObj->Self());
  int group;
  int groupRank;
  bool isOriginal = contract.isValid(self);
  if (isOriginal) {
    group = 0;
    groupRank = contract(self);
  } else {
    group = 1;
    /* MPI standard says they will be sorted by old rank,
       no need to spend time computing a good contiguous number */
    groupRank = 0;
  }
  MPI_Comm groupComm;
  MPI_Comm_split(MPI_COMM_WORLD, group, groupRank, &groupComm);
  CreateGroupCommResult result;
  result.isOriginal = isOriginal;
  result.group_pcu_obj = new pcu::PCU(groupComm);
  return result;
}


void getConfig(int argc, char** argv, pcu::PCU *PCUObj)
{
  if ( argc != 5 ) {
    if ( !PCUObj->Self() )
      printf("Usage: mpirun -n <outPartCount> %s"
             " <model> <inPartCount> <inMesh> <outMesh>\n"
             "Increase the part count of inMesh from inPartCount to outPartCount.\n"
             "Unlike the [z]split tool, outPartCount does not have to be an integer\n"
             "multiple of inPartCount.\n",
             argv[0]);
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }
  modelFile = argv[1];
  inputPartCount = atoi(argv[2]);
  meshFile = argv[3];
  outFile = argv[4];
  PCU_ALWAYS_ASSERT(inputPartCount <= PCUObj->Peers());
}

void balance(apf::Mesh2* m)
{
  apf::MeshTag* weights = m->createDoubleTag("zoltan_weight", 1);
  {
    apf::MeshIterator* it = m->begin(m->getDimension());
    apf::MeshEntity* e;
    double value = 1.0;
    while ((e = m->iterate(it)))
      m->setDoubleTag(e, weights, &value);
    m->end(it);
  }
  apf::Balancer* b = apf::makeZoltanBalancer(m, apf::RCB, apf::REPARTITION, false);
  b->balance(weights, 1.1);
  m->destroyTag(weights);
  delete b;
}

}


int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  {
  pcu::PCU expanded_pcu_obj = pcu::PCU(MPI_COMM_WORLD);
  lion_set_verbosity(1);
  gmi_register_mesh();
  getConfig(argc,argv,&expanded_pcu_obj);
  gmi_model* g = gmi_load(modelFile);
  apf::Mesh2* m = 0;
  CreateGroupCommResult result = createGroupComm(&expanded_pcu_obj);

  if (result.isOriginal)
    m = apf::loadMdsMesh(g, meshFile, result.group_pcu_obj);
  
  m = apf::expandMdsMesh(m, g, inputPartCount, &expanded_pcu_obj);
  balance(m);
  Parma_PrintPtnStats(m, "");
  m->writeNative(outFile);
  freeMesh(m);
  }
  MPI_Finalize();
}


