#include "LIIPBMod.h"

int main(int argc, char* argv[])
{
   ParUtil::Instance()->init(argc, argv);

   pGModel model = 0;
   pMesh mesh = MS_newMesh(model);
   
   PM_load(mesh, "geom.sms");

   int numParts = atoi(argv[1]);
   vector<pMesh> meshes;
   meshes.push_back(mesh);

   zoltanCB zlb;
   zlb.setLocalNumParts(numParts);
   M_loadbalance2(meshes,zlb);

   LIIPBMod liipbmod;
   liipbmod.run(meshes);

   ParUtil::Instance()->Finalize();


   return 1;
}
