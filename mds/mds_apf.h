/*
   Copyright 2014 Dan Ibanez

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
*/

#ifndef MDS_APF_H
#define MDS_APF_H

#ifdef __cplusplus
extern "C" {
#endif

#include "mds.h"
#include "mds_tag.h"
#include "mds_net.h"

struct mds_apf {
  struct mds mds;
  struct mds_tags tags;
  double (*point)[3];
  double (*param)[2];
  void** model[MDS_TYPES];
  void* user_model;
  void** parts[MDS_TYPES];
  struct mds_net remotes;
  struct mds_net matches;
};

struct mds_apf* mds_apf_create(void* model, int d, int cap[MDS_TYPES]);
void mds_apf_destroy(struct mds_apf* m);
double* mds_apf_point(struct mds_apf* m, mds_id e);
double* mds_apf_param(struct mds_apf* m, mds_id e);
void* mds_apf_model(struct mds_apf* m, mds_id e);
mds_id mds_apf_create_entity(
    struct mds_apf* m, int type, void* model, mds_id* from);
void mds_apf_destroy_entity(struct mds_apf* m, mds_id e);

void* mds_get_part(struct mds_apf* m, mds_id e);
void mds_set_part(struct mds_apf* m, mds_id e, void* p);

struct mds_tag* mds_number(struct mds_apf* m);
struct mds_apf* mds_reorder(struct mds_apf* m);

void* mds_find_model(struct mds_apf* m, int dim, int id);
int mds_model_dim(struct mds_apf* m, void* model);
int mds_model_id(struct mds_apf* m, void* model);

struct mds_apf* mds_read_smb(void* model, const char* pathname);
struct mds_apf* mds_write_smb(struct mds_apf* m, const char* pathname);

void mds_verify(struct mds_apf* m);
void mds_verify_residence(struct mds_apf* m, mds_id e);

extern int const mds_apf_double;
extern int const mds_apf_int;

#ifdef __cplusplus
}
#endif

#endif
