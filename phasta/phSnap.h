#ifndef PH_SNAP_H
#define PH_SNAP_H

#include <apf.h>
#include <apfMesh.h>
#include<stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

void pass_mesh_to_phasta(apf::Mesh* mesh);

void sim_get_pos_on_surf (double dx[], double dy[], double dz[], int numnp,
                          double px[], double py[], double pz[]);

void sim_is_in_closure (int e_dim, int e_tag, int t_dim, int t_tag, int& answer);

#ifdef __cplusplus
}
#endif

#endif
