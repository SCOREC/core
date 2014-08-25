#include <apf.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <gmi_mesh.h>
#include <parma.h>
#include <PCU.h>

apf::Migration* imbalanceVtx(apf::Mesh* m, int& goal) {
  apf::Migration* plan = new apf::Migration(m);
  if( PCU_Comm_Self() ) 
    return plan;
  int sent = 0;
  apf::MeshIterator* it = m->begin(0);
  apf::MeshEntity* e;
  while ( (e = m->iterate(it)) && sent < goal ) {
    if( ! m->isShared(e) ) 
      continue;
    apf::Copies rmts;
    m->getRemotes(e, rmts);
    apf::Adjacent elms;
    m->getAdjacent(e, m->getDimension(), elms);
    APF_ITERATE(apf::Adjacent, elms, elm) 
      if( ! plan->has(*elm) ) 
        plan->send(*elm, rmts.begin()->first);
    sent++;
  }
  m->end(it);
  goal -= sent;
  return plan;
}

int main(int argc, char** argv)
{
  assert(argc == 4);
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  PCU_Debug_Open();
  gmi_register_mesh();
  //load model and mesh
  apf::Mesh2* m = apf::loadMdsMesh(argv[1],argv[2]);
  int totTgt;
  int tgt = m->count(0)*.15*(!PCU_Comm_Self());
  do {
    apf::Migration* plan = imbalanceVtx(m, tgt);
    m->migrate(plan);
    totTgt = tgt;
    PCU_Max_Ints(&totTgt, 1);
    PCU_Debug_Print("tgt %d totTgt %d\n", tgt, totTgt);
  } while ( totTgt > 0 );
  m->writeNative(argv[3]);
  // destroy mds
  m->destroyNative();
  apf::destroyMesh(m);
  PCU_Comm_Free();
  MPI_Finalize();
}
