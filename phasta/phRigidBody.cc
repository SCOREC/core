#include <phastaChef.h>
#include <gmi_mesh.h>

#include "apfSIM.h"
#include "gmi_sim.h"
#include <SimUtil.h>
#include <SimModel.h>
#include <SimPartitionedMesh.h>
#include <SimMeshTools.h>

#include <pcu_util.h>
#include <cstdio>

#ifdef __cplusplus
extern"C"{
#endif

std::vector<ph::rigidBodyMotion> RBMs;

// update rbms
void core_update_rbms (double tx[], double ty[], double tz[],
                       double ax[], double ay[], double az[],
                       double px[], double py[], double pz[],
                       double ag[], double sc[],
                       int tags[],  int numRbm)
{
  RBMs.clear();
  for(int i = 0; i < numRbm; i++) {
    ph::rigidBodyMotion rbm = ph::rigidBodyMotion(tags[i], ag[i], sc[i]);
    rbm.set_trans(tx[i], ty[i], tz[i]);
    rbm.set_rotaxis(ax[i], ay[i], az[i]);
    rbm.set_rotpt(px[i], py[i], pz[i]);
    RBMs.push_back(rbm);
  }
}

// get rbms
void core_get_rbms (std::vector<ph::rigidBodyMotion> &rbms) {
  rbms = RBMs;
}

// get centroid of a model region
void core_get_centroid (int r_tag, double ct[]){
  apf::MeshSIM* apf_msim = dynamic_cast<apf::MeshSIM*>(m);
  pParMesh ppm = apf_msim->getMesh();
  pMesh pm = PM_mesh(ppm,0);

  gmi_model* gmiModel = apf_msim->getModel();
  pGModel model = gmi_export_sim(gmiModel);

  pGRegion modelRegion = (pGRegion) GM_entityByTag(model, 3, r_tag);
  pRegion meshRegion;
  double cent[3];
  double modelCent[3] = {0.0, 0.0, 0.0};
  double totalVol = 0.0;
  RIter rIter = M_classifiedRegionIter(pm, modelRegion);
  while((meshRegion=RIter_next(rIter))){
    EN_centroid(meshRegion, cent);
    double vol = R_volume(meshRegion);
    totalVol += vol;
    for (int i = 0; i < 3; i++)
      modelCent[i] += vol*cent[i];
  }
  RIter_delete(rIter);
  // divide by the total volume of the model region
  for (int i = 0; i < 3; i++)
    ct[i] = modelCent[i] / totalVol;
}

#ifdef __cplusplus
}
#endif

