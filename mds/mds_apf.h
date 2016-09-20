/****************************************************************************** 

  Copyright 2014 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/

#ifndef MDS_APF_H
#define MDS_APF_H

#include <gmi.h>

#ifdef __cplusplus
extern "C" {
#endif

#include "mds.h"
#include "mds_tag.h"
#include "mds_net.h"

struct pcu_file;

struct gmi_model;
struct gmi_ent;

struct mds_apf {
  struct mds mds;
  struct mds_tags tags;
  double (*point)[3];
  double (*param)[2];
  struct gmi_ent** model[MDS_TYPES];
  struct gmi_model* user_model;
  void** parts[MDS_TYPES];
  struct mds_net remotes;
//seol
  struct mds_net ghosts;
  struct mds_net matches;
};

struct mds_apf* mds_apf_create(struct gmi_model* model, int d,
    mds_id cap[MDS_TYPES]);
void mds_apf_destroy(struct mds_apf* m);
double* mds_apf_point(struct mds_apf* m, mds_id e);
double* mds_apf_param(struct mds_apf* m, mds_id e);
struct gmi_ent* mds_apf_model(struct mds_apf* m, mds_id e);
void mds_apf_set_model(struct mds_apf* m, mds_id e, struct gmi_ent* model);
mds_id mds_apf_create_entity(
    struct mds_apf* m, int type, struct gmi_ent* model, mds_id* from);
void mds_apf_destroy_entity(struct mds_apf* m, mds_id e);

void* mds_get_part(struct mds_apf* m, mds_id e);
void mds_set_part(struct mds_apf* m, mds_id e, void* p);

struct mds_tag* mds_number_verts_bfs(struct mds_apf* m);
struct mds_apf* mds_reorder(struct mds_apf* m, int ignore_peers,
    struct mds_tag* vert_numbers);

struct gmi_ent* mds_find_model(struct mds_apf* m, int dim, int id);
int mds_model_dim(struct mds_apf* m, struct gmi_ent* model);
int mds_model_id(struct mds_apf* m, struct gmi_ent* model);

struct mds_apf* mds_read_smb(struct gmi_model* model, const char* pathname,
    int ignore_peers, void* apf_mesh);
struct mds_apf* mds_write_smb(struct mds_apf* m, const char* pathname,
    int ignore_peers, void* apf_mesh);

void mds_verify(struct mds_apf* m);
void mds_verify_residence(struct mds_apf* m, mds_id e);

int mds_align_matches(struct mds_apf* m);
int mds_align_ghosts(struct mds_apf* m);
int mds_align_remotes(struct mds_apf* m);

void mds_derive_model(struct mds_apf* m);

extern int const mds_apf_double;
extern int const mds_apf_int;
extern int const mds_apf_long;

void mds_write_smb_meta(struct pcu_file* file, void* mesh_cpp);
void mds_read_smb_meta(struct pcu_file* file, struct mds_apf* mesh,
                       void* mesh_cpp);

#ifdef __cplusplus
}
#endif

#endif
