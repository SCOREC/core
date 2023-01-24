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

void core_get_surf_normal (int flag, int f_tag, int e_or_v_tag, double par1, double par2, double &n1, double &n2, double &n3) {

//input: flag, f_tag, e_tag, v_tag, par1, par2
//output: unit normals n1, n2, n3
  double fpar[2]; //parametric location to GF_normal
  double xyz[3];  //evaluated unit normal at the paramer location
  int f_dim = 2;
  int e_dim = 1;
  int v_dim = 0;
  n1 = 0;
  n2 = 0;
  n3 = 0;
  apf::ModelEntity* f  = m->findModelEntity(f_dim, f_tag);
//flag here is same value as m2gClsfcn(i,1), i.e., 0 = vertex, 1 = edge, 2 = face
  if(flag == 2) {//GFace
    fpar[0] = par1;
    fpar[1] = par2;
  }else if(flag == 1) { //GEdge
    apf::ModelEntity* e  = m->findModelEntity(e_dim, e_or_v_tag);
    if(GE_dirUsed(reinterpret_cast<pGEdge>(e), reinterpret_cast<pGFace>(f)) == 2){
      // std::cout << "core_get_surf_normal: Seam GEdge encountered, which is not allowed"<<std::endl;
      exit(EXIT_FAILURE);     
    }
    GF_edgeReparam(reinterpret_cast<pGFace>(f), reinterpret_cast<pGEdge>(e), par1, 0, fpar);	
  }else if(flag == 0) { //GVertex
    apf::ModelEntity* v  = m->findModelEntity(v_dim, e_or_v_tag);
    GF_vertexReparam(reinterpret_cast<pGFace>(f), reinterpret_cast<pGVertex>(v), fpar);
  }else{
    // std::cout << "core_get_surf_normal: Wrong flag is set. Allowed values are 1,2,3 only"<<std::endl;
    exit(EXIT_FAILURE);                  
  }

  GF_normal(reinterpret_cast<pGFace>(f), fpar, xyz);
  n1 = xyz[0];
  n2 = xyz[1];
  n3 = xyz[2];
}

#ifdef __cplusplus
}
#endif

