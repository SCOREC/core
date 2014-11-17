/****************************************************************************** 

  Copyright 2014 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/

#ifndef AGM_H
#define AGM_H

#ifdef __cplusplus
extern "C" {
#endif

struct agm;

enum agm_obj_type {
  AGM_ENTITY,
  AGM_USE,
  AGM_BOUNDARY,
  AGM_OBJ_TYPES
};

enum agm_ent_type {
  AGM_VERTEX,
  AGM_EDGE,
  AGM_FACE,
  AGM_REGION,
  AGM_ENT_TYPES
};

enum agm_use_type {
  AGM_VERTEX_USE,
  AGM_EDGE_USE,
  AGM_FACE_USE,
  AGM_USE_TYPES
};

enum agm_bdry_type {
  AGM_ENDPOINTS,
  AGM_LOOP,
  AGM_SHELL,
  AGM_BDRY_TYPES
};

struct agm_ent {
  enum agm_ent_type type;
  int id;
};

struct agm_use {
  enum agm_use_type type;
  int id;
};

struct agm_bdry {
  enum agm_bdry_type type;
  int id;
};

struct agm* agm_new(void);
void agm_free(struct agm* m);

struct agm_ent agm_add_ent(struct agm* m, enum agm_ent_type t);
struct agm_bdry agm_add_bdry(struct agm* m, struct agm_ent e);
struct agm_use agm_add_use(struct agm* m, struct agm_bdry b, struct agm_ent of);

void agm_reserve(struct agm* m, enum agm_ent_type t, int n);

int agm_ent_count(struct agm* m, enum agm_ent_type t);
int agm_use_count(struct agm* m, enum agm_use_type t);
int agm_bdry_count(struct agm* m, enum agm_bdry_type t);

struct agm_ent agm_first_ent(struct agm* m, enum agm_ent_type t);
struct agm_ent agm_next_ent(struct agm* m, struct agm_ent e);

int agm_ent_null(struct agm_ent e);
int agm_ent_eq(struct agm_ent a, struct agm_ent b);
int agm_use_null(struct agm_use u);
int agm_use_eq(struct agm_use a, struct agm_use b);
int agm_bdry_null(struct agm_bdry b);
int agm_bdry_eq(struct agm_bdry a, struct agm_bdry b);

struct agm_use agm_first_use_of(struct agm* m, struct agm_ent e);
struct agm_use agm_next_use_of(struct agm* m, struct agm_use u);
struct agm_use agm_first_use_by(struct agm* m, struct agm_bdry b);
struct agm_use agm_next_use_by(struct agm* m, struct agm_use u);
struct agm_ent agm_used(struct agm* m, struct agm_use u);
struct agm_bdry agm_user(struct agm* m, struct agm_use u);
struct agm_ent agm_bounds(struct agm* m, struct agm_bdry b);
struct agm_bdry agm_first_bdry_of(struct agm* m, struct agm_ent e);
struct agm_bdry agm_next_bdry_of(struct agm* m, struct agm_bdry b);

int agm_use_count_of(struct agm* m, struct agm_ent e);
int agm_use_count_by(struct agm* m, struct agm_bdry b);
int agm_bdry_count_of(struct agm* m, struct agm_ent e);
int agm_down_count(struct agm* m, struct agm_ent e);

struct agm_tag;

struct agm_tag* agm_new_tag(struct agm* m, int bytes);
void* agm_tag_at(struct agm_tag* t, enum agm_obj_type o,
    int subtype, int index);

enum agm_ent_type agm_type_from_dim(int dim);
int agm_dim_from_type(enum agm_ent_type t);

struct agm_use agm_find_use_by_bdry(struct agm* m, struct agm_ent of,
    struct agm_bdry by);
struct agm_use agm_find_use_by_ent(struct agm* m, struct agm_ent of,
    struct agm_ent by);
int agm_find_path(struct agm* m, struct agm_ent from, struct agm_ent to,
    struct agm_use path[4]);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif
