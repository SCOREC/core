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

apf::Mesh2* m;

ph::Input in;

void pass_info_to_phasta(apf::Mesh2* mesh, ph::Input& ctrl){
  m  = mesh;
  in = ctrl;
}

/* input dx,dy,dz are current coordinates of mesh in phastsa
   output px,py,pz are corrected coordinates of mesh */
void core_get_pos_on_surf (double dx[], double dy[], double dz[], int numnp,
                           int f[], double px[], double py[], double pz[]) {
  double newpar[2];
  double newpt[3];
  int counter = 0;
  apf::MeshEntity* v;
  apf::MeshIterator* it = m->begin(0);
  while ((v = m->iterate(it))) {
    apf::Vector3 p;
    m->getPoint(v, 0, p);
    pVertex meshVertex = reinterpret_cast<pVertex>(v);
    // only call this on face or edge mesh vertex
    if((V_whatInType(meshVertex)==2 || V_whatInType(meshVertex)==1)
       && (f[counter] == 1)) {
      const double disp[3] = {dx[counter]-p[0],dy[counter]-p[1],dz[counter]-p[2]};
      V_movedParamPoint(meshVertex,disp,newpar,newpt);
      px[counter] = newpt[0];
      py[counter] = newpt[1];
      pz[counter] = newpt[2];
    }
    else {
      px[counter] = dx[counter];
      py[counter] = dy[counter];
      pz[counter] = dz[counter];
    }
    counter++;
  }
  m->end(it);
  PCU_ALWAYS_ASSERT(counter==numnp);
}

void core_is_in_closure (int e_dim, int e_tag, int t_dim, int t_tag, int& answer) {
  answer = 0;
  apf::ModelEntity* e  = m->findModelEntity(e_dim, e_tag);
  apf::ModelEntity* et = m->findModelEntity(t_dim, t_tag);
  answer = m->isInClosureOf(e, et);
  PCU_ALWAYS_ASSERT(answer == 0 || answer == 1);
}

#ifdef __cplusplus
}
#endif

