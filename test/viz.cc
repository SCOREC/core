#include <gmi_mesh.h>
#include <apf.h>
#include <apfMesh.h>
#include <apfMesh2.h>
#include <apfMDS.h>
#include <PCU.h>
#include <parma.h>
#include "../viz/viz.h"

namespace {
  const char* modelFile = 0;
  const char* meshFile = 0;
  
  void freeMesh(apf::Mesh* m)
  {
    m->destroyNative();
    apf::destroyMesh(m);
  }

  void getConfig(int argc, char** argv)
  {
    assert(argc==3);
    modelFile = argv[1];
    meshFile = argv[2];
  }
  Color misColor(apf::Mesh* m) {
    int neighborDim = m->getDimension()-1;
    int misNumber = Parma_MisNumbering(m,neighborDim);
    Color colors[9] = 
      {RED,BLUE,GREEN,PURPLE,ORANGE,YELLOW,BROWN,PINK,GREY};
    return colors[misNumber%9];
  }
}
int main(int argc, char** argv)
{
  int provided;
  MPI_Init_thread(&argc,&argv,MPI_THREAD_MULTIPLE,&provided);
  assert(provided==MPI_THREAD_MULTIPLE);
  PCU_Comm_Init();
  gmi_register_mesh();
  getConfig(argc,argv);
  apf::Mesh2* m = apf::loadMdsMesh(modelFile,meshFile);
 
  
  Visualization v;
  
  char output[128];
  sprintf(output,"%d",PCU_Comm_Self());
  std::string part_num(output);

  apf::MeshIterator* itr;
  apf::MeshEntity* ent;
  v.watchMesh(m);
  v.breakpoint("The whole mesh");

  for (int i=0;i<3;i++) {
    v.watchDimension(m,2,misColor(m));
    v.watchDimension(m,0);
    itr = m->begin(1);
    while ((ent=m->iterate(itr))!=0) {
      if (m->isShared(ent))
        v.watchEntity(m,ent);

    }
    v.markPart(m,part_num);
    sprintf(output,"Testing MIS %d",i);
    v.breakpoint(std::string(output));
  }

  v.watchDimension(m,1,BYPART);
  v.breakpoint();

  Color c = misColor(m);
  itr=m->begin(2);
  int i=0;
  while((ent=m->iterate(itr))!=0) {
    if (i==0)
      v.watchEntity(m,ent,c);
    i++;
    i%=3;
  }

  v.breakpoint("Every third face");

  itr = m->begin(0);
  i=0;
  while((ent=m->iterate(itr))!=0) {
    if (i==0)
      v.watchEntity(m,ent,BYPART);
    i++;
    i%=2;
  }
  v.breakpoint();

  v.watchBoundary(m,1,BLACK);
  itr = m->begin(2);
  while ((ent=m->iterate(itr))!=0) {
    if (m->isShared(ent))
      v.watchEntity(m,ent);
  }
  
  v.showAxis();
  v.breakpoint("Part Boundaries");

  freeMesh(m);
  PCU_Comm_Free();
  MPI_Finalize();
}
