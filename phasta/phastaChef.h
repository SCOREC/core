#ifndef PHASTA_CHEF_H
#define PHASTA_CHEF_H

#include <apf.h>
#include <apfMesh.h>
#include <apfMesh2.h>
#include <phInput.h>
#include <stdlib.h>

namespace ph {
  void attachSIMSizeField(apf::Mesh2* m, apf::Field* sf_mag, apf::Field* sf_dir);
}

#ifdef __cplusplus
extern "C" {
#endif

extern apf::Mesh2* m;

extern ph::Input in;

void pass_info_to_phasta(apf::Mesh2* mesh, ph::Input ctrl);

void core_get_pos_on_surf (double dx[], double dy[], double dz[], int numnp,
                          double px[], double py[], double pz[]);

void core_is_in_closure (int e_dim, int e_tag, int t_dim, int t_tag, int& answer);

void core_measure_mesh (double x1[], double x2[], double x3[], int numnp,
                        double& minvq, double& minfq);

#ifdef __cplusplus
}
#endif

#endif
