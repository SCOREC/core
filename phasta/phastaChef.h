#ifndef PHASTA_CHEF_H
#define PHASTA_CHEF_H

#include <apf.h>
#include <apfMesh.h>
#include <apfMesh2.h>
#include <phInput.h>
#include <stdlib.h>
#include <cstring>
#include <vector>

namespace ph {
  void attachSIMSizeField(apf::Mesh2* m, apf::Field* sf_mag, apf::Field* sf_dir);
  void attachSIMSizeField(apf::Mesh2* m, apf::Field* sf_mag);

  struct rigidBodyMotion {
    rigidBodyMotion(int t = 0, double r = 0.0, double s = 0.0)
    {
      tag = t;
      memset(trans, 0.0, sizeof trans);
      memset(rotaxis, 0.0, sizeof rotaxis);
      memset(rotpt, 0.0, sizeof rotpt);
      rotang = r;
      scale = s;
    }
    int tag;
    double trans[3];
    double rotaxis[3];
    double rotpt[3];
    double rotang;
    double scale;
    void set_trans(double x, double y, double z){
      trans[0] = x; trans[1] = y; trans[2] = z;
    }
    void set_rotaxis(double x, double y, double z){
      rotaxis[0] = x; rotaxis[1] = y; rotaxis[2] = z;
    }
    void set_rotpt(double x, double y, double z){
      rotpt[0] = x; rotpt[1] = y; rotpt[2] = z;
    }
  };

} // end namespace

#ifdef __cplusplus
extern "C" {
#endif

extern apf::Mesh2* m;

extern ph::Input in;

void pass_info_to_phasta(apf::Mesh2* mesh, ph::Input& ctrl);

void core_get_pos_on_surf (double dx[], double dy[], double dz[], int numnp,
                           int f[], double px[], double py[], double pz[]);

void core_is_in_closure (int e_dim, int e_tag, int t_dim, int t_tag, int& answer);

void core_get_surf_normal (int flag, int f_tag, int e_or_v_tag, double par1, double par2, double& n1, double& n2, double& n3);

void core_ShockVelSmoother(double x1[], double x2[], double x3[], int numnp, 
                         double V_i[], double V_o[]);

void core_measure_mesh (double x1[], double x2[], double x3[], int numnp,
                        double& minvq, double& minfq);

void core_update_rbms (double tx[], double ty[], double tz[],
                         double ax[], double ay[], double az[],
                         double px[], double py[], double pz[],
                         double ag[], double sc[],
                         int tags[],  int numRbm);

void core_get_rbms (std::vector<ph::rigidBodyMotion> &rbms);

void core_get_centroid (int r_tag, double ct[]);

#ifdef __cplusplus
}
#endif

#endif
