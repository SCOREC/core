/****************************************************************************** 

  Copyright 2014 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/

#ifndef MDS_H
#define MDS_H

#include "mds_config.h"

enum {
  MDS_VERTEX,
  MDS_EDGE,
  MDS_TRIANGLE,
  MDS_QUADRILATERAL,
  MDS_WEDGE,
  MDS_PYRAMID,
  MDS_TETRAHEDRON,
  MDS_HEXAHEDRON,
  MDS_TYPES
};

typedef MDS_ID_TYPE mds_id;

#define MDS_NONE -1

struct mds {
  int d;
  mds_id n[MDS_TYPES];
  mds_id cap[MDS_TYPES];
  mds_id end[MDS_TYPES];
  int mrm[4][4];
  mds_id* down[4][MDS_TYPES];
  mds_id* up[4][MDS_TYPES];
  mds_id* first_up[4][MDS_TYPES];
  mds_id* free[MDS_TYPES];
  mds_id first_free[MDS_TYPES];
};

struct mds_set {
  int n;
  mds_id e[MDS_SET_MAX];
};

extern int const mds_dim[MDS_TYPES];
extern int const mds_degree[MDS_TYPES][4];
extern int const* mds_types[MDS_TYPES][4];

void mds_create(struct mds* m, int d, mds_id cap[MDS_TYPES]);
void mds_destroy(struct mds* m);
mds_id mds_create_entity(struct mds* m, int type, mds_id *from);
void mds_destroy_entity(struct mds* m, mds_id e);
mds_id mds_find_entity(struct mds* m, int type, mds_id *from);
int mds_type(mds_id e);
mds_id mds_index(mds_id e);
mds_id mds_identify(int type, mds_id idx);
void mds_get_adjacent(struct mds* m, mds_id e, int dim, struct mds_set* s);
mds_id mds_begin(struct mds* m, int dim);
mds_id mds_next(struct mds* m, mds_id);

void mds_add_adjacency(struct mds* m, int from_dim, int to_dim);
void mds_remove_adjacency(struct mds* m, int from_dim, int to_dim);

int mds_has_up(struct mds* m, mds_id e);

void mds_change_dimension(struct mds* m, int d);

#endif
